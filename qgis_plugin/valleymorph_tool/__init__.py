def classFactory(iface):
    from .plugin import ValleyMorphPlugin
    return ValleyMorphPlugin(iface)
