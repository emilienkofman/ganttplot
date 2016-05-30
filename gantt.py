#!/usr/bin/env python
# gantt.py ganttfile | gnuplot -p

import itertools
import argparse
from ConfigParser import ConfigParser
import sys, os, stat
from collections import OrderedDict
from pint import UnitRegistry
import json
from itertools import combinations
from pprint import pprint
import logging as log
log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)

ur = UnitRegistry()
ur.default_format = '~'


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
    def __init__(self, resource, start, stop, label, precedences=[], color="0xFFFFFF", border=1, bordercolor="0xFFFFFF", alpha=1, cumulative=0):
        self.resource = resource
        self.start = start
        self.stop = stop
        self.label = label
        self.precedences = precedences
        self.color = color
        self.border = border
        self.bordercolor = bordercolor
        self.alpha = alpha
        self.cumulative = cumulative

    def __str__(self):
        return self.label
    def __repr__(self):
        return str(self)

    def overlap(self, other):
        return (self.cumulative == other.cumulative and self.resource == other.resource and not (self.stop<=other.start or other.stop <= self.start))

    def rectangle(self):
        resmargin=0.1
        rescenter = self.resource.index
        resbottom = rescenter-0.5+resmargin
        restop = rescenter+0.5-resmargin
        resheight = restop-resbottom

        recstep = resheight/self.resource.cumulative
        reccenter = resbottom+recstep/2.+recstep*self.cumulative
        recbottom = reccenter-recstep/2.
        rectop = reccenter+recstep/2.

        bottomleft = (self.start, recbottom)
        topright = (self.stop, rectop)

        return Rectangle(bottomleft, topright, self.color, self.border, self.bordercolor, self.alpha)

    def arrows(self):
        arrows = []
        rfrom = self.rectangle()
        for target in self.precedences:
            rto = target.rectangle()
            if rfrom.under(rto):
                xyfrom = rfrom.topright
                xyto = rto.bottomleft
            elif rto.under(rfrom):
                xyfrom = rfrom.bottomright
                xyto = rto.topleft
            else:
                xyfrom = rfrom.topright[0], rfrom.center[1]
                xyto = rto.topleft[0], rto.center[1]
            arrows.append(Arrow(xyfrom, xyto))
        return arrows



class Rectangle(object):
    """
    Container for rectangle information.
    """
    rectangleHeight = 0.8  #: Height of a rectangle in units.
    def __init__(self, bottomleft, topright, fillcolor, border, bordercolor, alpha=1):
        self.bottomleft = bottomleft
        self.bottomright = topright[0], bottomleft[1]
        self.topright = topright
        self.topleft = bottomleft[0], topright[1]
        self.center = ((bottomleft[0]+topright[0])/2, (bottomleft[1]+topright[1])/2)
        self.fillcolor = fillcolor
        self.fillstyle = 'solid 0.8'
        self.linewidth = border
        self.bordercolor = bordercolor
        self.alpha = alpha

    def under(self, other):
        if self.center[1] < other.center[1]:
            return True
        return False


class Arrow(object):
    """
    Container for precedence arrow information
    """
    def __init__(self, xyfrom, xyto):
        self.xyfrom = xyfrom
        self.xyto = xyto


class Resource(object):
    """
    Container for resources information
    """
    def __init__(self, label, cumulative=1, index=0):
        self.label = label
        self.index = index
        self.cumulative = cumulative

def remove_overlaps(activities):
    for a in activities:
        # find every other activity that overlaps a
        overlaps = [a1 for a1 in activities  if a.overlap(a1)]
        log.info(overlaps)
        if len(overlaps)>1:
            resource = overlaps[0].resource
            if resource.cumulative < len(overlaps):
                resource.cumulative = len(overlaps)
            for i,o in enumerate(overlaps):
                o.cumulative = i

    for a0, a1 in itertools.combinations(activities, 2):
        log.info((a0, a1, a0.overlap(a1), a0.cumulative, a1.cumulative))

    return activities


def merge_activities_same_label(activities):
    change = True
    while change:
        change=False
        merged = {}
        for a0, a1 in combinations(activities, 2):
            if a0.resource==a1.resource and a0.label==a1.label and a0.color==a1.color:
                if a0.stop==a1.start:
                    a01 = Activity(a0.resource, a0.start, a1.stop, a0.label, a0.precedences, a0.color)
                    merged[a01] = (a0, a1)
                    break
                elif a1.stop==a0.start:
                    a10 = Activity(a0.resource, a1.start, a0.stop, a0.label, a0.precedences, a0.color)
                    merged[a10] = (a1, a0)
                    break
        for a, (a0, a1) in merged.iteritems():
            if a0 in activities: activities.remove(a0)
            if a1 in activities: activities.remove(a1)
            activities.append(a)
            change=True
    return activities


def parse_json_gantt(ganttlines, options):
    """
    Load the resource/task json file
    """
    activities = {}
    precedences = {}
    schedule = json.loads("\n".join(ganttlines))

    for activityID,data in schedule.iteritems():
        resource = data['map']
        if options.pinttime:
            start = ur(data['start'])
            stop = ur(data['stop'])
        else:
            start = float(data['start'])
            stop = float(data['stop'])
        label = data['label']
        precedences[activityID] = data['precedences']
        activities[activityID] = \
            Activity(resource, start, stop, label,
                     color=data.get("color","0xFFFFFF"),
                     border=int(data.get("border", "1")),
                     alpha=float(data.get("alpha", "1")),
                     bordercolor=data.get("bordercolor", "0x000000"))

    for activityID, activity in activities.iteritems():
        activity.precedences = [activities[aid] for aid in precedences[activityID]]


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
    log.info(resources)
    labels = [a.label for a in activities]

    # Sort such that resources and tasks appear in alphabetical order
    if alphasort:
        resources.sort()
        labels.sort()

    # Resources are read from top (y=max) to bottom (y=1)
    resources.reverse()
    resources = {r:Resource(r,index=i) for i,r in enumerate(resources)}

    for a in activities:
        a.resource = resources[a.resource]

    return labels, resources.values(), activities


def generate_plotdata(activities, resources, tasks, options):
    """
    Generate Gnuplot lines.
    """
    xmin = 0
    if not options.xmax:
        xmax = max(act.stop for act in activities)
    else:
        xmax = options.xmax
    ymin = -0.5
    ymax = len(resources)-0.5
    xlabel = options.xlabel
    ylabel = ''
    title = options.plottitle
    resourcenames = [r.label for r in resources]
    ytics = ''.join(['(',
                     ', '.join(('"%s" %d' % (r.label, r.index))
                                for r in resources),
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
                       'set xlabel "%s" font ",%d"' % (xlabel, options.fontsize),
                       'set tics font ",12"',
                       'set ylabel "%s" font ",%d"' % (ylabel, options.fontsize),
                       'set title "%s"' % title,
                       'set ytics %s' % ytics,
                       'set key %s' % key_position,
                       'set grid %s' % grid_tics,
                       'unset colorbox',
                       'set termopt enhanced']

    # Generate gnuplot rectangle objects
    plot_rectangles = []
    for n, a in enumerate(activities):
        r = a.rectangle()
        rectangle = ['set object %d rectangle' % (n+1),
                                 'from %f, %0.1f' % r.bottomleft,
                                 'to %f, %0.1f' % r.topright,
                                 'fillcolor rgb ' + r.fillcolor,
                                 'linecolor rgb '+r.bordercolor, # TODO: This won't work
                                 'linewidth '+str(r.linewidth)]

        style="fillstyle transparent solid "+str(r.alpha)
        if r.linewidth==0:
            style+=" noborder"

        rectangle.append(style)
        plot_rectangles.append(' '.join(rectangle))

    for a in activities:
        log.info(a.arrows())

    arrows = []
    for a in activities:
        arrows+=a.arrows()
    log.info(arrows)

    plot_arrows = (' '.join(['set arrow',
                                 'from %f, %0.1f' % a.xyfrom,
                                 'to %f, %0.1f' % a.xyto,
                                 'lw 2'])
                    for a in arrows)

    # Generate labels inside the rectangles
    rectangles = [a.rectangle() for a in activities]
    tcut = [(t[:10] + '..') if len(t) > 10 else t for t in tasks]
    plot_labels = (' '.join(['set label "%s"' %  t,
                             'at %f,' % r.center[0],
                             '%f' % r.center[1],
                             'rotate by +90',
                             'center font ",%d"' % options.fontsize])
                  for r, t in zip(rectangles, tcut))

    # Generate gnuplot lines
    plot_lines = ['plot ' +
                  ', \\\n\t'.join(' '.join(['-1',
                                      'title "%s"' % t,
                                      'with lines',
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



def pint_to_float(activities, options):
    """
    Convert python pint start/stop dates to floats,
    and get an xlabel such that it is human readable.
    """
    end = max(a.stop for a in activities)
    unit = end.to_compact().units
    options.xlabel = "time ("+str(unit)+")"
    for a in activities:
        a.start, a.stop = a.start.to(unit).magnitude, a.stop.to(unit).magnitude


def compute(options, ganttlines):
    activities = parse_json_gantt(ganttlines, options)

    if options.mergelabels:
        activities = merge_activities_same_label(activities)

    if options.pinttime:
        pint_to_float(activities, options)

    tasks, resources, activities = make_unique_tasks_resources(options.alphasort, activities)
    activities = remove_overlaps(activities)
    generators = generate_plotdata(activities, resources, tasks, options)

    write_data(generators, options)


parser = argparse.ArgumentParser(description='Transform a list of intervals associated with resources into a gantt diagram.')
parser.add_argument("-o", "--output", type=str, help='output filename', 
            default='', dest='outputfile')
parser.add_argument("-a", "--alphasort", help='', action="store_true", default=False)
parser.add_argument("--mergelabels", help='merge adjacent tasks with same labels', action="store_true", default=False)
parser.add_argument("-p", "--pinttime", help='Use python pint unit system for start/stop times', action="store_true",
            default=False)
parser.add_argument("-t", "--title", type=str, help='Title', 
            default='', dest='plottitle')
parser.add_argument("-l", "--xlabel", type=str, help='x label', 
            default='time', dest='xlabel')
parser.add_argument("-x", "--xmax", type=float, help='Fixed plot time')
parser.add_argument("-f", "--fontsize", type=int, help='', default=12)
parser.add_argument('--legend', dest='legend', action='store_true')
parser.add_argument('--no-legend', dest='legend', action='store_false')
parser.set_defaults(legend=False)

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


