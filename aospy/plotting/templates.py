class MultTemplates(object):
    def __init__(self):
        pass

# Annual cycle of monthly mean data from one model and run as maps.
acmm = MultTemplates()
acmm.an_row = 4
acmm.n_col = 3
acmm.n_panel = 12
acmm.n_data = 1
acmm.plot_type = 'map'
acmm.intvls = range(1,13)
acmm.ax_titles = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
