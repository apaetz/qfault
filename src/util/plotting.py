'''
Utility functions for plotting with matplotlib.

Created on 2010-10-24

@author: adam
'''
from counting.countParallel import iterParallel
import logging
import matplotlib
matplotlib.use('PDF')  # Save plots as PDF files.
from matplotlib import pyplot as plt

logger = logging.getLogger('plotting')

#fontMgr = matplotlib.font_manager.FontManager(size='x-small')
#fontProps = matplotlib.font_manager.FontProperties()


#'-'	solid line style
#'--'	dashed line style
#'-.'	dash-dot line style
#':'	dotted line style
#'.'	point marker
#','	pixel marker
#'o'	circle marker
#'v'	triangle_down marker
#'^'	triangle_up marker
#'<'	triangle_left marker
#'>'	triangle_right marker
#'1'	tri_down marker
#'2'	tri_up marker
#'3'	tri_left marker
#'4'	tri_right marker
#'s'	square marker
#'p'	pentagon marker
#'*'	star marker
#'h'	hexagon1 marker
#'H'	hexagon2 marker
#'+'	plus marker
#'x'	x marker
#'D'	diamond marker
#'d'	thin_diamond marker
#'|'	vline marker
#'_'	hline marker
markers = ['o', '^', 's', '*', 'D', 'p', 'h', 'x']




#'b'	blue
#'g'	green
#'r'	red
#'c'	cyan
#'m'	magenta
#'y'	yellow
#'k'	black
#'w'	white
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']





def evalExprList(expr, X):
	'''
	Evaluates expression expr for all values in the list X.
	'''
	return [r.get() for r in iterParallel(X, expr)]
	#return [expr(x) for x in X]

def evalExpr(val, expr):
	'''
	Evaluates expression expr at value val.  Useful for parallel evaluation.
	'''
	return expr(val)

def lineStyle(i):
	marker = markers[i % len(markers)]
	color = getColor(i)
	return str(color) + str(marker) + '-'

def getColor(i, colors=colors):
	return colors[i % len(colors)]

class Series(object):
	
	def __init__(self, X, Y, name):
		self.X = X
		self.Y = Y
		self.name = name


def plotInit():
	plt.hold(False)
	plt.figure()
	plt.hold(True)
	
def plotShow(filename=None):
	plt.show()
	if filename != None:
		plt.savefig(filename)
	
	plt.hold(False)

def barSeries(seriesList, filename=None, legendLoc='upper right', width=0.75, colors=colors):
	X = set()
	for i,series in enumerate(seriesList):
		X |= set(series.X)
		plt.bar(series.X, series.Y, width=width, color=getColor(i,colors), label=series.name)
		
	w2 = width/2.
	X = sorted(X)
	plt.xticks([x + w2 for x in X], [str(x) for x in X])
	
	plt.legend(loc=legendLoc)
	plotShow(filename)
	

def plotPolyList(polyList, xMin, xMax, 
			filename=None, 
			labelList=None, 
			xLabel='', 
			yLabel='', yscale='linear',  ylim=None,
			legendLoc='upper right', legendCols=1, legendPoints=2, legendFrame=True,
			numPoints=100):
	
	
	plotInit()
	
	plt.ylabel(yLabel)
	plt.xlabel(xLabel)
	
	xStep = (xMax-xMin)/(numPoints-1)
	X = [xMin + i*xStep for i in range(numPoints)]
	
	handles = []
	for e, poly in enumerate(polyList):
		Y = evalExprList(poly, X)
		#Y = [poly(x) for x in X]
		
		ls = lineStyle(e)
		if None != labelList:
			print labelList[e], Y
			handles.append(plt.plot(X,Y, ls, label=labelList[e]))
		else:
			print e+1, Y			
			handles.append(plt.plot(X,Y, ls))
	

	# Place legend at the right side of the axes viewport
	#plt.figlegend(handles, labelList, loc='upper right')
	plt.legend(loc=legendLoc, ncol=legendCols, numpoints=legendPoints,frameon=legendFrame)
	plt.yscale(yscale)
	if None != ylim:
		plt.ylim(ylim)

	plotShow(filename)
	
def plotList(X, yList, yErrList=None, filename=None, labelList=None, xLabel='', yLabel='', legendLoc='upper right',
			 xscale='linear', yscale='linear', xlim=None, ylim=None,
			 custom_command=None):
	
	#fontProps.set_size('small')
	plt.hold(False)
	plt.figure()
	plt.hold(True)
	#plt.axes([0.2, 0.3, 0.65, 0.65])
	
	plt.ylabel(yLabel)
	plt.xlabel(xLabel)
	
	if None == yErrList:
		yErrList = [None] * len(yList)
		
	if None == labelList:
		labelList = [None] * len(yList)
	
	handles = []
	for e, Y in enumerate(yList):
		
		ls = lineStyle(e)
		print labelList[e], Y, yErrList[e]
		handles.append(plt.errorbar(X,Y, yerr=yErrList[e], fmt=ls, label=labelList[e]))		



	# Place legend at the right side of the axes viewport
	#plt.figlegend(handles, labelList, loc='upper right')
	plt.legend(loc=legendLoc)
	
	plt.xscale(xscale)
	plt.yscale(yscale)
	
	if None != xlim:
		plt.xlim(xlim)
	if None != ylim:
		plt.ylim(ylim)
		
	if None != custom_command:
		custom_command(plt)
		
#	yaxis = plt.gca().yaxis
#	interval = yaxis.get_view_interval()
#	print 'interval=', interval
#	yaxis.set_view_interval(0, interval[1])
#	print 'interval=', yaxis.get_view_interval()

	plt.show()
	if filename != None:
		plt.savefig(filename)
	
	plt.hold(False)



def plot(poly, xMin, xMax, numPoints=100, label=None):
	from counting.countParallel import iterParallel
	
	dx = (xMax - xMin) / numPoints
	X = [xMin + dx*i for i in range(numPoints)]
	results = iterParallel(X, evalExpr, [poly])
	Y = [r.get() for r in results]
	
	plt.plot(X,Y, label=label)
	plt.legend()
	plt.show()