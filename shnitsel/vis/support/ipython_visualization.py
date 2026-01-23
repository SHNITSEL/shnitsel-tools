import abc

from matplotlib.figure import Figure, SubFigure

from shnitsel.vis.support.visualizeable import Visualizable


class IPythonVisualizable:
    def __init__(self, source: Figure | SubFigure | Visualizable):
        """
        # Use HTML to visualize multiple outputs in IPython output:

        from IPython.display import display, HTML
        def make_html(folder, image):
            return '<img src="{}" style="display:inline;margin:1px"/>'
                    .format(os.path.join(folder, image))

        display(HTML(''.join(make_html(f, x)) for x in files))

        # More in-depth options using IPython widgets:
        from IPython.display import display
        from ipywidgets import widgets, HBox

        imageA = widgets.Image(value=open('path/to/image/a.jpg', 'rb').read())
        imageB = widgets.Image(value=open('path/to/image/b.jpg', 'rb').read())

        hbox = HBox([imageA, imageB])
        display(hbox)

        # optional: you can show more hboxes, boxes will be arranged vertically
        # display(anotherHbox) # 2nd
        # display(yetAnotherHbox) # 3rd
        # display(andAnotherHbox) # 4th
        """
