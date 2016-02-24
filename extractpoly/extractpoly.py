#!/usr/bin/env python 
'''
Copyright (C) 2026 Bart Pieters, b.pieters@fz-juelich.de

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
'''

'''
B. Pieters 18.01.2016
This file is based on Aaron Spike's flatten Bezier and code snippets
It is mdified in that I now output the result to a text file rather than 
applying it. Furthermore I convert units to a selected unit.
'''
import inkex, cubicsuperpath, simplepath, cspsubdiv, sys
import os

class MyEffect(inkex.Effect):
    def __init__(self):
        inkex.Effect.__init__(self)
        self.OptionParser.add_option("-f", "--flatness",
                        action="store", type="float", 
                        dest="flat", default=10.0,
                        help="Minimum flatness of the subdivided curves")
        self.OptionParser.add_option("--filepath", action="store",
                                        type="string", dest="filepath",
                                        default=None, help="Filename and path to stora the polygon data")
        self.OptionParser.add_option("--unit", action="store",
                                        type="string", dest="unit",
                                        default=None, help="Length unit for the polygon data")
    def effect(self):
                
	# get page size:
        root_doc = self.document.getroot()
	viewbox = root_doc.attrib['viewBox']
	# width = root_doc.attrib['width']
	# height = root_doc.attrib['height']
	splitbox = viewbox.split(' ')
	fconv = float(self.unittouu(self.options.unit))
	x1 = float(splitbox[0])
	y1 = float(splitbox[1])
	# x2 = float(splitbox[2])
	y2 = float(splitbox[3])
	# w = float(self.unittouu(width))
	# h = float(self.unittouu(height))
	
	# Flatten path(s) and format polygon as a string
	outputString = ""
        for id, node in self.selected.iteritems():
            if node.tag == inkex.addNS('path','svg'):
                d = node.get('d')
                p = cubicsuperpath.parsePath(d)
                cspsubdiv.cspsubdiv(p, self.options.flat)
                for sp in p:
                    first = True
                    for csp in sp:
                        cmd = 'L'
                        if first:
                            cmd = 'M'
                        first = False
                        outputString += str((csp[1][0]-x1)/fconv) + "\t" + str((y2-csp[1][1]+y1)/fconv) + "\n"
                    outputString += "\n"
    
        # get the filename to write to
        path = self.options.filepath
        if (path != None ):
            if (not os.path.isabs(path)):
                if os.name == 'nt':
                    path = os.path.join(os.environ['USERPROFILE'],path)
                else:
                    path = os.path.join(os.path.expanduser("~"),path)
                f=open(path,'w')
		f.write(outputString)
		f.close
                inkex.errormsg('polygon extracted to: %s' % path)
	else:	
		outputString = ("# No filename was provided, using stderr\n") + outputString
                sys.stderr.write(outputString)

if __name__ == '__main__':
    e = MyEffect()
    e.affect()


# vim: expandtab shiftwidth=4 tabstop=8 softtabstop=4 fileencoding=utf-8 textwidth=99
