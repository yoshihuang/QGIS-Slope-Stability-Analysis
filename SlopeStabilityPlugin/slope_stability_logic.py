## -*- coding: utf-8 -*-
"""
QGIS Plugin Logic for Slope Stability Analysis (Bishop Method)
Integrates General Scan and Event-Based Analysis.
"""

import json
import numpy as np
import traceback
from qgis.PyQt.QtCore import QThread, pyqtSignal
from shapely.wkt import loads as wkt_loads
from shapely.geometry import LineString, Point, Polygon, mapping
from shapely.ops import unary_union

class SlopeStabilityWorker(QThread):
    # Define signals: progress(percent, message), finished(file_path), error(error_message)
    progress = pyqtSignal(int, str)
    finished = pyqtSignal(str)
    error = pyqtSignal(str)

    def __init__(self, geo_layer, water_layer, kh, mode, output_path):
        super().__init__()
        self.geo_layer = geo_layer
        self.water_layer = water_layer
        self.kh = kh
        self.mode = mode  # "general", "event", "both"
        self.output_path = output_path
        
        # Analysis results container
        self.results = []
        
        # Stratigraphy and geometry cache
        self.layers = []
        self.total_slope_poly = None
        self.bounds = None
        self.slope_direction = 1
        self.water_geom = None

    def run(self):
        """Thread entry point"""
        try:
            self.progress.emit(5, "Reading and parsing QGIS layers...")
            
            # 1. Parse geological model (Polygon)
            # Convert QGIS Features to Python dicts and Shapely objects to avoid thread conflicts
            if not self.geo_layer:
                raise ValueError("Geological layer not specified!")

            valid_features_count = 0
            for feat in self.geo_layer.getFeatures():
                geom = feat.geometry()
                if geom.isEmpty(): continue
                
                # Convert geometry
                shapely_poly = wkt_loads(geom.asWkt())
                
                # Extract attributes (attempt automatic casting)
                try:
                    props = {
                        'cohesion': float(feat['cohesion']),
                        'phi': float(feat['phi']),
                        'unit_weight': float(feat['unit_weight'])
                    }
                    self.layers.append({'poly': shapely_poly, 'props': props})
                    valid_features_count += 1
                except KeyError as e:
                    raise KeyError(f"Missing layer fields: {e}. Please ensure the layer contains 'cohesion', 'phi', and 'unit_weight' fields.")
                except ValueError:
                    raise ValueError(f"Invalid attribute value (non-numeric). Feature ID: {feat.id()}")

            if valid_features_count == 0:
                raise ValueError("No valid geological data found in the layer.")

            # Merge geometries to get the outer boundary
            self.total_slope_poly = unary_union([l['poly'] for l in self.layers])
            self.bounds = self.total_slope_poly.bounds # (minx, miny, maxx, maxy)
            self.slope_direction = self._detect_slope_direction()

            # 2. Parse water table (LineString) - Used only in Event mode if selected
            if self.water_layer:
                water_lines = []
                for f in self.water_layer.getFeatures():
                    geom = f.geometry()
                    if not geom.isEmpty():
                        water_lines.append(wkt_loads(geom.asWkt()))
                if water_lines:
                    self.water_geom = unary_union(water_lines)

            # 3. Start analysis
            total_tasks = 0
            if self.mode in ["general", "both"]: total_tasks += 1
            if self.mode in ["event", "both"]: total_tasks += 1
            
            current_task = 0

            # --- Mode A: General Scan ---
            if self.mode in ["general", "both"]:
                current_task += 1
                progress_start = 10
                progress_end = 50 if self.mode == "both" else 90
                self.progress.emit(progress_start, "Executing general matrix scan...")
                
                self.run_general_analysis(progress_range=(progress_start, progress_end))

            # --- Mode B: Event Analysis ---
            if self.mode in ["event", "both"]:
                current_task += 1
                progress_start = 50 if self.mode == "both" else 10
                progress_end = 90
                self.progress.emit(progress_start, f"Executing single event analysis (Kh={self.kh})...")
                
                self.run_event_analysis(progress_range=(progress_start, progress_end))

            # 4. Output results
            self.progress.emit(95, "Writing GeoJSON...")
            self._save_geojson()
            self.finished.emit(self.output_path)

        except Exception as e:
            # Catch all errors and report to main window
            error_msg = f"Analysis failed:\n{str(e)}\n\n{traceback.format_exc()}"
            self.error.emit(error_msg)

    def _detect_slope_direction(self):
        """Detect slope direction: 1 (Left High, Right Low), -1 (Right High, Left Low)"""
        mid_x = (self.bounds[0] + self.bounds[2]) / 2
        left_box = Polygon([(self.bounds[0], self.bounds[1]), (mid_x, self.bounds[1]), 
                            (mid_x, self.bounds[3]), (self.bounds[0], self.bounds[3])])
        left_area = self.total_slope_poly.intersection(left_box).area
        total_area = self.total_slope_poly.area
        return 1 if left_area > (total_area / 2) else -1

    def _get_pore_pressure(self, x, y_arc, water_method, water_param):
        """
        Calculate pore water pressure u (kPa)
        water_method: 'depth' or 'geometry'
        """
        wl_height = -9999.0
        
        if water_method == 'depth':
            # Simple method: Assume water table parallels surface at depth 'water_param'
            # Simplified here as deducting from model top
            # (Advanced version could find surface y at x then subtract depth)
            mid_x = (self.bounds[0] + self.bounds[2]) / 2
            slope_tan = (self.bounds[3] - self.bounds[1]) / (self.bounds[2] - self.bounds[0]) * -self.slope_direction
            # Base altitude set to slope top minus depth
            base_alt = self.bounds[3] - water_param
            wl_height = base_alt + slope_tan * (x - mid_x)
            
        elif water_method == 'geometry' and self.water_geom is not None:
            # Geometry method: Find water table elevation at x
            v_line = LineString([(x, self.bounds[1] - 100), (x, self.bounds[3] + 100)])
            isect = self.water_geom.intersection(v_line)
            
            if not isect.is_empty:
                if isect.geom_type == 'Point':
                    wl_height = isect.y
                elif isect.geom_type == 'MultiPoint':
                    # If multiple intersections (e.g., folded strata), take the highest point
                    wl_height = max([p.y for p in isect.geoms])
                elif isect.geom_type in ['GeometryCollection']:
                     for g in isect.geoms:
                         if g.geom_type == 'Point':
                             wl_height = max(wl_height, g.y)
        
        # Calculate u = gamma_w * hw
        if wl_height > y_arc:
            return 9.81 * (wl_height - y_arc)
        else:
            return 0.0

    def calculate_fs_stable(self, center, radius, kh, water_method, water_param):
        """Bishop Simplified Method Core Engine"""
        # 1. Geometric slicing
        angles = np.linspace(0, 2 * np.pi, 300)
        circle = LineString([(center[0] + radius * np.cos(a), center[1] + radius * np.sin(a)) for a in angles])
        clipped = circle.intersection(self.total_slope_poly)
        
        if clipped.is_empty: return None
        
        # Handle multiple intersections, keep the longest segment
        target_line = clipped
        if clipped.geom_type in ['MultiLineString', 'GeometryCollection']:
            max_len = 0
            for g in clipped.geoms:
                if g.geom_type == 'LineString' and g.length > max_len:
                    max_len = g.length
                    target_line = g
            if max_len == 0: return None
        elif clipped.geom_type != 'LineString':
            return None

        coords = list(target_line.coords)
        if len(coords) < 12: return None

        # 2. Outcrop check (relaxed to 5m tolerance)
        p_start, p_end = Point(coords[0]), Point(coords[-1])
        boundary = self.total_slope_poly.boundary
        if boundary.distance(p_start) > 5.0 or boundary.distance(p_end) > 5.0:
            return None
        # Bottom protection
        if p_start.y <= self.bounds[1] + 1.0 or p_end.y <= self.bounds[1] + 1.0:
            return None

        # 3. Slice integration
        res_sum, dri_sum = 0.0, 0.0
        slice_data = []

        for i in range(len(coords)-1):
            p1, p2 = coords[i], coords[i+1]
            x_mid, dx = (p1[0] + p2[0])/2, abs(p2[0] - p1[0])
            if dx < 0.05: continue
            
            y_arc = (p1[1] + p2[1])/2
            
            # Find surface height (y_top)
            v_line = LineString([(x_mid, self.bounds[1]-10), (x_mid, self.bounds[3]+500)])
            y_isect = self.total_slope_poly.intersection(v_line)
            if y_isect.is_empty: continue
            y_top = y_isect.bounds[3]
            
            # Angle calculation (lock to sliding direction)
            dy = p1[1] - p2[1] if self.slope_direction == 1 else p2[1] - p1[1]
            dX = p2[0] - p1[0] if self.slope_direction == 1 else p1[0] - p2[0]
            alpha = np.arctan2(dy, dX)
            
            if abs(np.degrees(alpha)) > 82: continue # Vertical filter

            # Get layer properties at this point
            mid_pt = Point(x_mid, y_arc)
            c, phi, gamma = 5.0, 20.0, 19.0
            
            # Simple Point in Polygon test
            found_layer = False
            for layer in self.layers:
                if layer['poly'].contains(mid_pt):
                    c = layer['props']['cohesion']
                    phi = layer['props']['phi']
                    gamma = layer['props']['unit_weight']
                    found_layer = True
                    break
            
            # If point is on boundary, contains might fail, use distance
            if not found_layer:
                for layer in self.layers:
                    if layer['poly'].distance(mid_pt) < 0.1:
                        c = layer['props']['cohesion']
                        phi = layer['props']['phi']
                        gamma = layer['props']['unit_weight']
                        break

            # Slice weight
            h_slice = max(0, y_top - y_arc)
            weight = h_slice * dx * gamma
            
            # Pore water pressure
            u = self._get_pore_pressure(x_mid, y_arc, water_method, water_param)

            # Effective normal force (no tension)
            ca, sa = np.cos(alpha), np.sin(alpha)
            # Prevent cos(alpha) from approaching 0
            safe_ca = ca if abs(ca) > 0.01 else 0.01 * (1 if ca >=0 else -1)
            
            n_eff = max(0, weight * ca - u * dx / abs(safe_ca))
            
            # Fellenius Term (as initial value)
            # res = c * (dx/safe_ca) + n_eff * np.tan(np.radians(phi))
            # dri = weight * sa + kh * weight * ca
            # res_sum += max(0, res)
            # dri_sum += dri
            
            slice_data.append({
                'w': weight, 'alpha': alpha, 'c': c, 'phi': np.radians(phi), 
                'u_dx': u * dx, 'dx': dx
            })

        # If no valid slices
        if not slice_data: return None

        # 4. Bishop iteration
        # Initial FS guess
        fs = 1.5 
        
        for _ in range(15):
            r_iter, d_iter = 0.0, 0.0
            for s in slice_data:
                # m_alpha = cos(a) + sin(a)tan(phi)/FS
                m_alpha = np.cos(s['alpha']) + (np.sin(s['alpha']) * np.tan(s['phi'])) / fs
                m_alpha = max(0.2, m_alpha) # Numerical protection
                
                term1 = s['c'] * s['dx'] + (s['w'] - s['u_dx']) * np.tan(s['phi'])
                r_iter += term1 / m_alpha
                
                # Driving Force
                d_iter += (s['w'] * np.sin(s['alpha']) + kh * s['w'] * np.cos(s['alpha']))
            
            if d_iter <= 1.0: return None # Driving force too small
            
            new_fs = r_iter / d_iter
            if abs(new_fs - fs) < 0.002: 
                fs = new_fs
                break
            fs = new_fs

        # Final return
        return float(fs), target_line

    def _execute_search(self, kh, water_method, water_param, bounds_cfg, steps, scenario_name, progress_cb):
        """Shared search logic"""
        cx_ext, cy_range, r_factor = bounds_cfg
        mw = self.bounds[2] - self.bounds[0]
        mh = self.bounds[3] - self.bounds[1]
        
        x_grid = np.linspace(self.bounds[0] - cx_ext[0]*mw, self.bounds[2] + cx_ext[1]*mw, steps[0])
        y_grid = np.linspace(self.bounds[3] + cy_range[0], self.bounds[3] + cy_range[1], steps[1])
        r_grid = np.linspace(mh * r_factor[0], mw * r_factor[1], steps[2])
        
        total_steps = len(x_grid) * len(y_grid) * len(r_grid)
        current_step = 0
        
        # Report progress every N steps to avoid GUI freeze
        report_interval = max(1, total_steps // 10) 

        for cx in x_grid:
            for cy in y_grid:
                for r in r_grid:
                    current_step += 1
                    if current_step % report_interval == 0:
                        pct = int((current_step / total_steps) * 100)
                        # This is just an internal micro-progress; macro progress bar is controlled externally
                        pass

                    if cy - r > self.bounds[3]: continue
                    
                    res = self.calculate_fs_stable((cx, cy), r, kh, water_method, water_param)
                    if res:
                        fs, geom = res
                        # Save only reasonable FS
                        if 0.1 < fs < 10.0:
                            self.results.append({
                                "type": "Feature",
                                "properties": {
                                    "Scenario": scenario_name,
                                    "Kh": float(kh),
                                    "FS": round(fs, 3),
                                    "Radius": round(r, 1)
                                },
                                "geometry": mapping(geom)
                            })

    def run_general_analysis(self, progress_range):
        """General Scan: Kh=[0, 0.1, 0.2], Depth=[0, 15, 30]"""
        kh_list = [0.0, 0.1, 0.2]
        depths = [0, 15, 30]
        
        total_scenarios = len(kh_list) * len(depths)
        completed = 0
        p_min, p_max = progress_range
        
        # Set sparser grid for speed
        # cx_ext, cy_range, r_factor
        bounds_cfg = ((0.5, 0.5), (20, 250), (0.5, 2.0)) 
        steps = (10, 8, 8) 

        for kh in kh_list:
            for d in depths:
                scenario = f"General_Kh{kh}_D{d}m"
                self._execute_search(kh, 'depth', d, bounds_cfg, steps, scenario, None)
                
                completed += 1
                current_p = p_min + int((p_max - p_min) * (completed / total_scenarios))
                self.progress.emit(current_p, f"General Scan: {scenario}")

    def run_event_analysis(self, progress_range):
        """Event Analysis: Specified Kh and Water Table"""
        p_min, p_max = progress_range
        
        # Set high-density grid
        bounds_cfg = ((0.5, 0.5), (10, 300), (0.5, 2.5))
        steps = (20, 15, 15)
        
        water_method = 'geometry' if self.water_geom else 'depth'
        water_param = self.water_geom if self.water_geom else 0.0 # Default depth 0 if no water layer
        
        scenario = f"Event_Kh{self.kh}"
        if water_method == 'depth': scenario += "_NoWaterLayer"
        else: scenario += "_WithWaterLayer"
        
        self.progress.emit(p_min, f"Starting high-precision analysis: {scenario}...")
        self._execute_search(self.kh, water_method, water_param, bounds_cfg, steps, scenario, None)
        self.progress.emit(p_max, "Event analysis completed")

    def _save_geojson(self):
        if not self.results:
            raise RuntimeError("Analysis completed, but no valid slip surfaces found.\nPossible reasons: insufficient search range or incorrect layer properties.")
            
        geojson_obj = {
            "type": "FeatureCollection",
            "name": "Slope_Stability_Results",
            "features": self.results
        }
        
        with open(self.output_path, 'w', encoding='utf-8') as f:
            json.dump(geojson_obj, f, indent=2)