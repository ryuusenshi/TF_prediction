import pandas
from math import pi

from bokeh.io import show
from bokeh.models import (
    ColumnDataSource,
    HoverTool,
    LinearColorMapper,
    BasicTicker,
    PrintfTickFormatter,
    ColorBar,
)
from bokeh.plotting import figure
#from bokeh.sampledata.unemployment1948 import data


df = pandas.read_csv('matrix_genotypes.txt', delimiter='\t')
df = df.drop(['_probePairKey'], axis=1)
df.index = df.index.map(str)
markers = list(df.index)
conditions = list(df.columns)

df = pandas.DataFrame(df.stack(), columns=['SNP']).reset_index()
df = df.rename(columns={'level_0': 'Marker', 'level_1': 'Condition'})


# this is the colormap from the original NYTimes plot
colors = ["#75968f", "#dfccce", "#550b1d"]
mapper = LinearColorMapper(palette=colors, low=0, high=2)

source = ColumnDataSource(df)

TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

p = figure(title="Genotype data",
           x_range=markers, y_range=conditions,
           x_axis_location="above", plot_width=900, plot_height=400,
           tools=TOOLS, toolbar_location='below')

p.grid.grid_line_color = None
p.axis.axis_line_color = None
p.axis.major_tick_line_color = None
p.axis.major_label_text_font_size = "5pt"
p.axis.major_label_standoff = 0
p.xaxis.major_label_orientation = pi / 3
p.xaxis.axis_label = 'Markers'
p.yaxis.axis_label = 'Conditions'
#p.xaxis.ticker = BasicTicker(desired_num_ticks=100)
#p.yaxis.ticker = BasicTicker(desired_num_ticks=100)

p.rect(x="Marker", y="Condition", width=1, height=1,
       source=source,
       fill_color={'field': 'SNP', 'transform': mapper},
       line_color=None)

color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="20pt",
                     ticker=BasicTicker(desired_num_ticks=2),
                     formatter=PrintfTickFormatter(format="%d"),
                     label_standoff=3, border_line_color=None, location=(0, 0))
p.add_layout(color_bar, 'right')

show(p)      # show the plot
