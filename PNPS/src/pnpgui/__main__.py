if __name__ == '__main__':
    import wx
    from . import pnpgui
    import logging
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    app = wx.App()
    pnpgui.pApp = app
    top = pnpgui.PNPFrame(None, None)
    top.Show()
    app.MainLoop()
