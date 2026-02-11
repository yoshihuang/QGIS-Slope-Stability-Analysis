def classFactory(iface):
    from .slope_stability_plugin import SlopeStabilityPlugin
    return SlopeStabilityPlugin(iface)