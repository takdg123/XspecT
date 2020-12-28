def AddGRB(ra,dec):
	Name='<source name="GRB" type="PointSource">\n'
	spec='\t<spectrum type="PowerLaw2">\n'
	spec+='\t\t<parameter free="1" max="1000.0" min="1e-05" name="Integral" scale="1e-06" value="1.0"/>\n'
	spec+='\t\t<parameter free="1" max="-1.0" min="-5.0" name="Index" scale="1.0" value="-2.0"/>\n'
	spec+='\t\t<parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="100.0"/>\n'# %epiv
	spec+='\t\t<parameter free="0" max="200000.0" min="20.0" name="UpperLimit" scale="1.0" value="10000.0"/>\n' 
	spec+='\t</spectrum>\n'
	skydir='\t<spatialModel type="SkyDirFunction">\n'
	skydir+='\t\t<parameter free="0" max="360.0" min="-360.0" name="RA" scale="1.0" value="%s"/>\n' %ra
	skydir+='\t\t<parameter free="0" max="90" min="-90" name="DEC" scale="1.0" value="%s"/>\n' %dec
	skydir+='\t</spatialModel>\n'
	skydir+='</source>'
	(src,)=(Name+spec+skydir,)
	return src