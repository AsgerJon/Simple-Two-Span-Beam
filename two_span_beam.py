# Two Span Beam Master File
# Import this as is to the notebook
# Invoke by writing: 
# *** CODE ***
# from two_span_beam import * 
# show
# *** END CODE ***
# in a cell and run the cell

from cont_beam_two_spans02 import *
from beam_drawing import *
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets

N = 1001

resetknap = widgets.ToggleButton(
    value=False,
    description='Klik her for at begynde',
    disabled=False,
    button_style='', # 'success', 'info', 'warning', 'danger' or ''
    tooltip='Description',
    icon='calculator', # (FontAwesome names without the `fa-` prefix)
    layout = widgets.Layout(
        width = '600px',
        height = '50px',
        ),
)

last_interval = widgets.IntRangeSlider(
    value=[0, 100],
    min=0,
    max=100,
    step=1,
    description='Last Interval som Procent af Bjælkens Længde:',
    disabled=False,
    continuous_update=False,
    orientation='horizontal',
    readout=True,
    readout_format='d',
    layout = widgets.Layout(
        width = '600px',
        height = '50px',
        ),
    style = {'description_width': 'initial',
             'handle_color': '#010000',
            },
)

last = widgets.FloatText(
    value = 1,
    description = 'Last i [N/mm]',
    disabled = False,
    layout = widgets.Layout(
        width = '600px',
        height = '25px',
        ),
    style = {'description_width': 'initial',
             'handle_color': '#010000',
            },
)

length = widgets.FloatText(
    value = 2000,
    description = 'Bjælkens Længde i [mm]:',
    disabled = False,
    layout = widgets.Layout(
        width = '600px',
        height = '25px',
        ),
    style = {'description_width': 'initial',
             'handle_color': '#010000',
            },
)

def last_funk(onoff, x, L, q):
    a0 = min([x[0],50])/100*L
    b0 = max([x[1],50])/100*L
    if onoff: draw_beam(a0, L-b0, L, q, N=1001)

show = widgets.interactive(last_funk, 
	onoff=resetknap,
	x=last_interval, 
	L=length, 
	q=last,
)