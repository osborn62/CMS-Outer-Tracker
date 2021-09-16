
class mpa_configurations():
	def __init__(self):
		self.nrowsraw = 16
		self.ncolsraw = 120
		self.npixsraw = 1920
		self.rowOneraw = 1 #not that row 0 does NOT exist
		self.colOneraw = 2 #not that column 0 and 1 still exist
		self.rowsraw = range(1,17)
		self.colsraw = range(0,120)
		self.colsraweff = range(2,120)

		self.nrowsnom = 16
		self.ncolsnom = 118
		self.npixsnom = 1888
		self.rowsnom = range(0,16) #shift by 1
		self.colsnom = range(0,118)#drop col 0,1, and shift by 2
		
		self.delaypulsebefore = 500
		self.delaypulseafter  = 100
		self.delayresetafter  = 200

	  
	def pixelidraw(self,row,column):
		return (row-1)*120+column  
	def colfrompixraw(self,pixelid):
		return pixelid%120
	def rowfrompixraw(self,pixelid):
		return pixelid/120 + 1

	def pixelidnom(self,row,column):
		return row*118+column	
	def colfrompixnom(self,pixelid):
		return pixelid%118
	def rowfrompixnom(self,pixelid):
		return pixelid/118

	def getrawrow(self,row):
		return row+1
	def getnomrow(self,row):
		return row-1
	def getrawcol(self,col):
		return col+2
	def getnomcol(self,col):
		return col-2
	def getrawpix(self,pix):
		row = self.getrawrow(self.rowfrompixnom(pix))
		col = self.getrawcol(self.colfronpixnom(pix))
		return pixelidraw(row,col)
	def getnompix(self,pix):
		row = self.getnomrow(self.rowfrompixraw(pix))
		col = self.getnomcol(self.colfronpixraw(pix))
		return pixelidnom(row,col)
	def convertRawToNomPixmap(self,data):
		if isinstance(data, list) and len(data)<self.npixsraw:
			print("Input is not a raw pixel map for sure, too short list - return original list")
			return data
		nomdata = []
		for p in range(0,self.npixsraw):
			#rows counting are identical, but columns with index 0,1 do not exist, drop them
			if self.colfrompixraw(p)==0 or self.colfrompixraw(p)==1:
				continue
			nomdata.append(data[p])
		return nomdata
	def convertNomToRawPixmap(self,data):
		if len(data)<self.npixsnom:
			print("Input is not a nom pixel map for sure, too short list - return original list")
			return data
		rawdata = []
		for p in range(0,self.npixsraw):
			#rows counting are identical, but columns with index 0,1 do not exist in nom, add them
			if self.colfrompixraw(p)==0 or self.colfrompixraw(p)==1:
				rawdata.append(0)
			else:
				nomp = self.getnompix(p)
				rawdata.append(data[nomp])
				return nomdata
	def getNomAndRawRowsCols(self,rows,cols,israw=False):
		wasempty = False
		if len(rows)==0:
			rows = self.rowsnom
			wasempty = True
		if len(cols)==0:
			cols = self.colsnom
			wasempty = True
		rawrows, rawcols = [],[]
		rawrows = [r for r in rows] 
		rawcols = [c for c in cols]
		if wasempty:
			rawrows = self.rowsraw
			rawcols = self.colsraw
		elif israw==False:
			for index in range(len(rawrows)): rawrows[index] = self.getrawrow(rows[index])
			for index in range(len(rawcols)): rawcols[index] = self.getrawcol(cols[index])
		else:
			for index in range(len(rows)): rows[index] = self.getnomrow(rows[index])
			for index in range(len(cols)): cols[index] = self.getnomcol(cols[index])
		return rows,cols,rawrows,rawcols
	def getNomAndRawRowCol(self,row,col,israw=False):
		rawrow,rawcol = 0,0
		rawrow = row
		rawcol = col
		if israw==False:
			rawrow = self.conf.getrawrow(row)
			rawcol = self.conf.getrawcol(col)
		else:
			row = self.conf.getnomrow(row)
			col = self.conf.getnomcol(col)
		return row,col,rawrow,rawcol
	def getPercentage(self, pixels):
		return str(round(100*float(len(pixels))/float(self.npixsnom),2))+"%"


conf = mpa_configurations()
