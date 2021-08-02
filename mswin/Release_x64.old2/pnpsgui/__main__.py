if __name__ == '__main__':
    import wx
    from . import pnpsgui
    import logging
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    app = wx.App()
    pnpsgui.pApp = app
    top = pnpsgui.PNPFrame(None, None)
    top.Show()
    app.MainLoop()
