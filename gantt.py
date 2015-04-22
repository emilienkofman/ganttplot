#!/usr/bin/env python
#
# TODO:
# - Task colors:
#     - User-defined using config file.
#     - Automagically chosen from color space.
#     - Advanced algorithm (contact Hannes Pretorius).
# - Koos' specs:
#     - Resources and tasks sorted in read-in order (default)
#       or alphabetically (flag).
#     - Have proper gnuplot behavior on windows/x11, eps/pdf, latex terminals.
#     - Create and implement algorithm for critical path analysis.
# - Split generic stuff into a Gantt class, and specific stuff into the main.
#
# gantt.py ganttfile | gnuplot

import itertools
import argparse
from ConfigParser import ConfigParser
import sys, os, stat
from collections import OrderedDict
from pprint import pprint

rectangleHeight = 0.8  #: Height of a rectangle in units.

class Activity(object):
    """
    Container for activity information.

    @ivar resource: Resource name.
    @type resource: C{str}

    @ivar start: Start time of the activity.
    @type start: C{float}

    @ivar stop: End time of the activity.
    @type stop: C{float}

    @ivar task: Name of the task/activity being performed.
    @type task: C{str}
    """
    def __init__(self, resource, start, stop, task, precedences=[]):
        self.resource = resource
        self.start = start
        self.stop = stop
        self.task = task
        self.precedences = precedences
    def __str__(self):
        return self.task+"["+str(self.start)+"-"+str(self.stop)+"]"
    def __repr__(self):
        return str(self)

class Rectangle(object):
    """
    Container for rectangle information.
    """
    def __init__(self, bottomleft, topright, fillcolor):
        self.bottomleft = bottomleft
        self.topright = topright
        self.center = ((bottomleft[0]+topright[0])/2, (bottomleft[1]+topright[1])/2)
        self.fillcolor = fillcolor
        self.fillstyle = 'solid 0.8'
        self.linewidth = 2

class Arrow(object):
    """
    Container for precedence arrow information
    """
    def __init__(self, xyfrom, xyto):
        self.xyfrom = xyfrom
        self.xyto = xyto


class ColorBook(object):
    """
    Class managing colors.

    @ivar colors
    @ivar palette
    @ivar prefix
    """
    def __init__(self, colorfname, tasks):
        """
        Construct a ColorBook object.

        @param colorfname: Name of the color config file (if specified).
        @type  colorfname: C{str} or C{None}

        @param tasks: Existing task types.
        @type  tasks: C{list} of C{str}
        """
        if colorfname:
            values = self.load_config(colorfname, tasks)
        else:
            values = self.fixed(tasks)

        self.colors, self.palette, self.prefix = values


    def load_config(self, colorfname, tasks):
        """
        Read task colors from a configuration file.
        """
        palettedef = 'model RGB'
        colorprefix = 'rgb'

        # Read in task colors from configuration file
        config = ConfigParser()
        config.optionxform = str # makes option names case sensitive
        config.readfp(open(colorfname, 'r'))
        # Colors are RGB colornames
        colors = dict(config.items('Colors'))

        # Raise KeyError if no color is specified for a task
        nocolors = [t for t in tasks if not colors.has_key(t)]
        if nocolors:
            msg = 'Could not find task color for ' + ', '.join(nocolors)
            raise KeyError(msg)

        return colors, palettedef, colorprefix

    def fixed(self, tasks):
        """
        Pick colors from a pre-defined palette.
        """
        # Set task colors
        # SE colors
        # (see http://w3.wtb.tue.nl/nl/organisatie/systems_engineering/\
        #      info_for_se_students/how2make_a_poster/pictures/)
        # Decrease the 0.8 values for less transparent colors.
        se_palette = {"se_red":   (1.0, 0.8, 0.8),
                     "se_pink":   (1.0, 0.8, 1.0),
                     "se_violet": (0.8, 0.8, 1.0),
                     "se_blue":   (0.8, 1.0, 1.0),
                     "se_green":  (0.8, 1.0, 0.8),
                     "se_yellow": (1.0, 1.0, 0.8)}
        se_gradient = ["se_red", "se_pink", "se_violet",
                       "se_blue", "se_green", "se_yellow"]
        se_palettedef = '( ' + \
                        ', '.join(('%d ' % n +
                                   ' '.join((str(x) for x in se_palette[c]))
                                   for n, c in enumerate(se_gradient))) + \
                        ' )'

        palettedef = 'model RGB defined %s' % se_palettedef
        colorprefix = 'palette frac'
        # Colors are fractions from the palette defined
        if len(tasks)-1!=0:
            colors = dict((t, '%0.2f' % (float(n)/(len(tasks)-1)))
                       for n, t in enumerate(tasks))
        else:
            colors = {tasks[0]: '0'}

        return colors, palettedef, colorprefix

def make_rectangles(activities, resource_map, colors):
    """
    Construct a collection of L{Rectangle} for all activities.

    @param activities: Activities being performed.
    @type  activities: C{iterable} of L{Activity}

    @param resource_map: Indices of all resources.
    @type  resource_map: C{dict} of C{str} to C{int}

    @param colors: Colors for all tasks.
    @type  colors: C{dict} of C{str} to C{str}

    @return: Collection of rectangles to draw.
    @rtype:  C{list} of L{Rectangle}
    """
    rectangles = []
    for act in activities:
        ypos = resource_map[act.resource]
        bottomleft = (act.start, ypos - 0.5 * rectangleHeight)
        topright = (act.stop, ypos + 0.5 * rectangleHeight)
        fillcolor = colors[act.task]
        rectangles.append(Rectangle(bottomleft, topright, fillcolor))

    return rectangles


def make_arrows(activities, resource_map, colors):

    arrows = []
    for source in activities:
        sourcecenter = resource_map[source.resource]
        sourceright = source.stop
        sourcebottom = sourcecenter - 0.5 * rectangleHeight
        sourcetop = sourcecenter + 0.5 * rectangleHeight

        for target in source.precedences:
            targetcenter = resource_map[target.resource]
            targetleft = target.start
            targetbottom = targetcenter - 0.5 * rectangleHeight
            targettop = targetcenter + 0.5 * rectangleHeight
            if sourcecenter > targetcenter:
                xyfrom = (sourceright, sourcebottom)
                xyto = (targetleft, targettop)
            elif sourcecenter < targetcenter:
                xyfrom = (sourceright, sourcetop)
                xyto = (targetleft, targetbottom)
            elif sourcecenter == targetcenter:
                xyfrom = (sourceright, sourcetop)
                xyto = (targetleft, targettop)

            arrows.append(Arrow(xyfrom, xyto))

    return arrows


def parse_gantt(ganttlines):
    """
    Load the resource/task file.

    @param ganttfile: Name of the gantt file.
    @type  ganttfile: C{str}

    @return: Activities loaded from the file, collection of
             (resource, start, end, task) activities.
    @rtype:  C{list} of L{Activity}
    """
    activities = {}
    precedences = {}
    for line in ganttlines: 
        line = line.strip().split()
        if len(line) == 0:
            continue
        resource = line[0]
        start = float(line[1])
        stop = float(line[2])
        task = line[3]
        precedences[task] = line[4:]
        activities[task] = Activity(resource, start, stop, task)
    
    for task, a in activities.iteritems():
        a.precedences = [activities[task] for task in precedences[a.task]]

    return activities.values()


def make_unique_tasks_resources(alphasort, activities):
    """
    Construct collections of unique task names and resource names.

    @param alphasort: Sort resources and tasks alphabetically.
    @type  alphasort: C{bool}

    @param activities: Activities to draw.
    @type  activities: C{list} of L{Activity}

    @return: Collections of task-types and resources.
    @rtype:  C{list} of C{str}, C{list} of C{str}
    """
    # Create list with unique resources and tasks in activity order.
    resources = list(OrderedDict.fromkeys([a.resource for a in activities]))
    tasks = list(OrderedDict.fromkeys([a.task for a in activities]))

    # Sort such that resources and tasks appear in alphabetical order
    if alphasort:
        resources.sort()
        tasks.sort()

    # Resources are read from top (y=max) to bottom (y=1)
    resources.reverse()

    return tasks, resources


def generate_plotdata(activities, resources, tasks, rectangles, arrows, options,
                     resource_map, color_book):
    """
    Generate Gnuplot lines.
    """
    xmin = 0
    if not options.xmax:
        xmax = max(act.stop for act in activities)
    else:
        xmax = options.xmax
    ymin = 0 + (rectangleHeight / 2)
    ymax = len(resources) + 1 - (rectangleHeight / 2)
    xlabel = 'time'
    ylabel = ''
    title = options.plottitle
    ytics = ''.join(['(',
                     ', '.join(('"%s" %d' % item)
                                for item in resource_map.iteritems()),
                     ')'])
    # outside and 2 characters from the graph
    if options.legend:
        key_position = 'outside width +2'
    else:
        key_position = 'off'
    grid_tics = 'xtics'

    # Set plot dimensions
    plot_dimensions = ['set xrange [%f:%f]' % (xmin, xmax),
                       'set yrange [%f:%f]' % (ymin, ymax),
                       'set autoscale x', # extends x axis to next tic mark
                       'set xlabel "%s"' % xlabel,
                       'set ylabel "%s"' % ylabel,
                       'set title "%s"' % title,
                       'set ytics %s' % ytics,
                       'set key %s' % key_position,
                       'set grid %s' % grid_tics,
                       'set palette %s' % color_book.palette,
                       'unset colorbox',
                       'set termopt enhanced']

    # Generate gnuplot rectangle objects
    plot_rectangles = (' '.join(['set object %d rectangle' % n,
                                 'from %f, %0.1f' % r.bottomleft,
                                 'to %f, %0.1f' % r.topright,
                                 'fillcolor %s %s' % (color_book.prefix,
                                                      r.fillcolor),
                                 'fillstyle solid 0.8'])
                    for n, r in itertools.izip(itertools.count(1), rectangles))

    plot_arrows = (' '.join(['set arrow',
                                 'from %f, %0.1f' % a.xyfrom,
                                 'to %f, %0.1f' % a.xyto,
                                 'lw 2'])
                    for a in arrows)

    # Generate labels inside the rectangles
    plot_labels = (' '.join(['set label "%s"' %t,
                             'at %f,' % r.center[0],
                             '%f' % r.center[1],
                             'rotate by +90',
                             'center'])
                  for r, t in zip(rectangles, tasks))

    # Generate gnuplot lines
    plot_lines = ['plot ' +
                  ', \\\n\t'.join(' '.join(['-1',
                                      'title "%s"' % t,
                                      'with lines',
                                      'linecolor %s %s ' % (color_book.prefix,
                                                        color_book.colors[t]),
                                      'linewidth 6'])
                            for t in tasks)]

    return plot_dimensions, plot_rectangles, plot_labels, plot_arrows, plot_lines

def write_data(generators, options):
    """
    Write plot data out to file or screen.

    @param fname: Name of the output file, if specified.
    @type  fname: C{str}  (??)
    """
    if options.outputfile:
        g = open(options.outputfile, 'w')
        g.write('\n'.join(itertools.chain(*generators)))
        g.close()
    else:
        print '\n'.join(itertools.chain(*generators))

def fmt_opt(short, long, arg, text):
    if arg:
        return '-%s %s, --%s%s\t%s' % (short[:-1], arg, long, arg, text)
    else:
        return '-%s, --%s\t%s' % (short, long, text)


def compute(options, ganttlines):
    activities = parse_gantt(ganttlines)
    tasks, resources = make_unique_tasks_resources(options.alphasort,
                                                   activities)

    # Assign indices to resources
    resource_map = dict(itertools.izip(resources, itertools.count(1)))

    color_book = ColorBook(options.colorfile, tasks)
    rectangles = make_rectangles(activities, resource_map, color_book.colors)
    arrows     = make_arrows(activities, resource_map, color_book.colors)

    generators = generate_plotdata(activities, resources, tasks, rectangles, arrows,
                    options, resource_map, color_book)

    write_data(generators, options)


parser = argparse.ArgumentParser(description='Transform a list of intervals associated with resources into a gantt diagram.')
parser.add_argument("-o", "--output", type=str, help='output filename', 
            default='', dest='outputfile')
parser.add_argument("-c", "--color", type=str, help='colors filename', 
            default='', dest='colorfile')
parser.add_argument("-a", "--alphasort", help='', action="store_true", 
            default=False)
parser.add_argument("-t", "--title", type=str, help='Title', 
            default='', dest='plottitle')
parser.add_argument("-x", "--xmax", type=int, help='Fixed plot time')
parser.add_argument('--legend', dest='legend', action='store_true')
parser.add_argument('--no-legend', dest='legend', action='store_false')
parser.set_defaults(legend=True)


if __name__ == '__main__':
    mode = os.fstat(0).st_mode
    ganttlines = []
    if stat.S_ISFIFO(mode) or stat.S_ISREG(mode): 
        ganttlines = sys.stdin.readlines()
    else: # nothing is passed on stdin, expect a filename
        parser.add_argument("filename")
        options = parser.parse_args()
        with open(options.filename, 'r') as infile:
            ganttlines = infile.readlines()

    options = parser.parse_args()
    compute(options, ganttlines)


