
class Table(list):
    def __init__(self, initlist=[], header=[]):
        initlist = [list(header)]+list(initlist)
        super(Table, self).__init__(initlist)
    def _repr_html_(self):
        html = ["<table>"]
        for row in self:
            html.append("<tr>")
            for col in row:
                try:
                    if col > 1e-3 or (abs(col) < 1e-20):
                        html.append("<td><DIV align='right'>{:.3f}</DIV></td>".format(col))
                    else:
                        html.append("<td><DIV align='left'>{:.3e}</DIV></td>".format(col))
                except:
                    html.append("<td><DIV align='left'>{}</DIV></td>".format(col))
            html.append("</tr>")
        html.append("</table>")
        return ''.join(html) 

class FitTable(list):
    def __init__(self, initlist=[], parlist = [], LaTex=False):
        if LaTex:
            initlist = [["From", "To", "Model"]+parlist+["Stats", "dof", "BIC"]]+list(initlist)
        else:
            initlist = [["From", "To", "Model"]+parlist+["Stats", "dof", "BIC", "Itv #"]]+list(initlist)
        super(FitTable, self).__init__(initlist)
    def _repr_html_(self):
        html = ["<table align='center'>"]
        for row in self:
            html.append("<tr align='center'>")
            for col in row:
                html.append("<td><DIV align='left'>{}</DIV></td>".format(col))
            html.append("</tr>")
        html.append("</table>")
        return ''.join(html)
    def _repr_latex_(self):
        latex = ["\\begin{longtable*}"]
        latex.append("{"+" ".join((["c"]*len(self[0])))+"}\n")
        for row in self:
            latex.append(" & ".join(map(format, row)))
            latex.append("\\\\\hline \n")
        latex.append("\\end{longtable*}")
        return ''.join(latex)

class ModelTable(list):
    def __init__(self, initlist=[]):
        initlist = [["", "Model", "", "Model", "", "Model", "", "Model", "", "Model"]]+list(initlist)
        super(ModelTable, self).__init__(initlist)
    def _repr_html_(self):
        html = ["<table>"]
        for row in self:
            html.append("<tr>")
            for col in row:
                try:
                    int(col)
                    html.append("<td><DIV align='right'>{}:</DIV></td>".format(col))
                except:
                    html.append("<td><DIV align='left'>{}</DIV></td>".format(col))
            html.append("</tr>")
        html.append("</table>")
        return ''.join(html)   

class IntvTable(list):
    def __init__(self, initlist=[]):
        initlist = [["", "Time [s]", "", "Time [s]", "", "Time [s]", "", "Time [s]", "", "Time [s]"]]+list(initlist)
        super(IntvTable, self).__init__(initlist)
    def _repr_html_(self):
        html = ["<table>"]
        for row in self:
            html.append("<tr>")
            for col in row:
                try:
                    int(col)
                    html.append("<td><DIV align='right'>{}:</DIV></td>".format(col))
                except:
                    html.append("<td><DIV align='left'>{}</DIV></td>".format(col))
            html.append("</tr>")
        html.append("</table>")
        return ''.join(html)   

