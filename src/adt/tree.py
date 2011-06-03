'''
Created on 2011-04-29

@author: adam
'''

class TreeNode(object):
	
	def __init__(self, name, values, parent=None, children=[]):
		self._name = name
		self._values = values
		self._parent = parent
		self._children = {}
		for c in children:
			self.addChild(c)
		
	def name(self):
		return self._name
	
	def val(self, key):
		return self._values[key]

	def __getitem__(self, key):
		return self.val(key)
	
	def __setitem__(self, key, val):
		self._values[key] = val
	
	def parent(self):
		return self._parent
	
	def setParent(self, parent):
		oldParent = self.parent()
		self._parent = None
		if None != oldParent:
			oldParent.removeChild(self)
		self._parent = parent
	
	def child(self, name):
		return self._children[name]
	
	def children(self):
		return [self._children[name] for name in sorted(self._children.keys())]
	
	def addChild(self, child):
		self._children[child.name()] = child
		child.setParent(self)
		
	def removeChild(self, child):
		self._children.remove(child.name())
		child.setParent(None)
		
	def __str__(self):
		string = ''
		
		if None != self.parent():
			for _, val in sorted(self.parent()._values.items()):
				string += str(val) + '.'
			string += '<-' 
		
		string += self.name() + '='
		
#		for _, val in sorted(self._values.items()):
#			string += str(val) + '.'
#					
#		for c in sorted(self.children()):
#			string += str(c._values) + ':'

		string += self.subtreeStr()
		
		return string
	
	def subtreeStr(self):
		string = ''
		for _, val in sorted(self._values.items()):
			string += str(val) + '.'
			
		if 0 != len(string):
			string = string[:-1]
			
		children = self.children()
		if 0 != len(children):
			string += '->'
			
		for c in children:
			string += c.subtreeStr() + ':'
			
		if 0 != len(children):
			string = string[:-1]
			
		return string
