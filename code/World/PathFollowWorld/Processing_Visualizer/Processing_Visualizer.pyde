# this file will display the output generated by berry world 
# run in visualization mode.

# provide path to HarvestWorldData file in FileName
# left and right arrows alter fps
# left mouse pauses and restarts

### replace path ###
fileName = 'c:/Users/cliff//ALIFE_EXTENDED/data2/pathVisualization.txt'
#fileName = 'c:/Users/cliff//WORK/MABE_DEV_3_26_2021_FRAG/PathFollow_TEST_2_extra_out/pathVisualization.txt'
#fileName = 'c:/Users/cliff//fragWork2/Jan_31_path/C00_161_pathVisualization.txt'
#fileName = 'c:/Users/cliff//fragWork2/Jan_31_path/C05_170_pathVisualization.txt'
#fileName = 'c:/Users/cliff//fragWork2/Jan_31_path/C05_159_pathVisualization.txt'
#fileName = 'c:/Users/cliff//fragWork2/Jan_31_path/pathVisualization.txt'


showSignals = True   # dispay turn signals on each map
windowSizeX = 800  # change the height of the window
windowSizeY = 900  # change the width of the window

#small map size
#windowSizeX = 500  # change the height of the window
#windowSizeY = 600  # change the width of the window
gridSize = 25       # change scale

fps = 4 # Frames Per Second, can also be adjusted with left and right arrows

############################################################
##
##  do not edit anything below!
##
############################################################

# global vars
fileHandle = open(fileName, 'r+')

directions = 8
worldTime = 0
mapTime = 0

mode = 1

def mousePressed():
    global mode
    if mode == 1:
        noLoop()
        mode = 0
    elif mode == 0:
        loop()
        mode = 1

def keyPressed():
    global fps
    if key == CODED:
        if (keyCode == RIGHT):
            fps = fps * 2
            frameRate(fps)
        if (keyCode == LEFT):
            fps = fps / 2;
            if (fps<1):
                fps = 1;
            frameRate(fps);

def readNextLineFromFile():
    global fileHandle
    line = fileHandle.readline().strip()
    while (line == ""):
        line = fileHandle.readline().strip()
    splitLine = split(line,',')
    return splitLine

           
def setup():
    global fileHandle
    global fps
    global directions
    global mode
    global gridSize
    global windowSizeX
    global windowSizeY
        
    background(0)
    stroke(0)
    frameRate(fps)
    size(windowSizeX, windowSizeY)
    smooth(30)
    textSize(gridSize)
    mousePressed()
  
def draw():
    global fps
    global mode
    global gridSize
    global showSignals
    
    global worldTime
    global mapTime

    
    clear()
    splitLine = readNextLineFromFile()
    if splitLine[0] == "new": # first pass
        mapTime = 1
        splitLine = readNextLineFromFile() # pass new
    if splitLine[0] == "start": # first time on this map
        worldTime = 1
        splitLine = readNextLineFromFile() # pass start
        
    dir = int(splitLine[0])
    splitLine = readNextLineFromFile() #read score
    fill(255)
    
    text("map   worldTime          score",0,gridSize)
    text("  "+str(mapTime)+"              "+str(worldTime)+  "            " + splitLine[0],0,gridSize*2)
    splitLine = readNextLineFromFile() #read x and y size
    splitLine = readNextLineFromFile() #read x and y size

    splitLine = readNextLineFromFile() #read turn right signal
    rs = str(splitLine[0])
    splitLine = readNextLineFromFile() #read turn left signal
    ls = str(splitLine[0])
    if showSignals:
        text("Left Signal: "+ls+"  Right Signal: " + rs,0,gridSize*3)

    splitLine = readNextLineFromFile() #read first line of map
    
    yOffset = 4 # start drawing map with a y offset to leave room for text
    # draw the map (until we find a start (new image) or new (new map))
    while splitLine[0] != "start" and splitLine[0] != "new":
        xOffset = 0
        while xOffset < len(splitLine[0]):
            ch = splitLine[0][xOffset]
            if ch == '*':
                fill(255)
                rect(xOffset * gridSize, yOffset * gridSize, gridSize, gridSize)
                fill(255,0,0)
                if dir == 0:
                    rect(xOffset * gridSize + gridSize * .3, yOffset * gridSize + gridSize * 0, gridSize * .4, gridSize * .4)
                if dir == 1:
                    rect(xOffset * gridSize + gridSize * .6, yOffset * gridSize + gridSize * 0, gridSize * .4, gridSize * .4)
                if dir == 2:
                    rect(xOffset * gridSize + gridSize * .6, yOffset * gridSize + gridSize * .3, gridSize * .4, gridSize * .4)
                if dir == 3:
                    rect(xOffset * gridSize + gridSize * .6, yOffset * gridSize + gridSize * .6, gridSize * .4, gridSize * .4)
                if dir == 4:
                    rect(xOffset * gridSize + gridSize * .3, yOffset * gridSize + gridSize * .6, gridSize * .4, gridSize * .4)
                if dir == 5:
                    rect(xOffset * gridSize + gridSize * .0, yOffset * gridSize + gridSize * .6, gridSize * .4, gridSize * .4)
                if dir == 6:
                    rect(xOffset * gridSize + gridSize * .0, yOffset * gridSize + gridSize * .3, gridSize * .4, gridSize * .4)
                if dir == 7:
                    rect(xOffset * gridSize + gridSize * .0, yOffset * gridSize + gridSize * 0, gridSize * .4, gridSize * .4)

            elif ch == '1':
                fill(128,255,128)
                rect(xOffset * gridSize, yOffset * gridSize, gridSize, gridSize)
            #elif ch == 's':
            #    newValue = 0;
            #    xOffset+=1
            #    while(splitLine[xOffset] != 's'):
            #        newValue *= 10
            #        newValue += int(splitLine[j])
            #    xOffset+=1
                    
            elif ch == '4':
                fill(128,128,255)
                rect(xOffset * gridSize, yOffset * gridSize, gridSize, gridSize)
                fill(0)
                text('E',xOffset * gridSize + gridSize * .3, yOffset * gridSize + gridSize * .9)
            elif ch == '2':
                fill(210)
                ellipse((xOffset+.5) * gridSize, (yOffset+.5) * gridSize, gridSize, gridSize)
                fill(0)
                text(rs,xOffset * gridSize + gridSize * .2, yOffset * gridSize + gridSize * .9)
            elif ch == '3':
                fill(210)
                ellipse((xOffset+.5) * gridSize, (yOffset+.5) * gridSize, gridSize, gridSize)
                fill(0)
                text(ls,xOffset * gridSize + gridSize * .2, yOffset * gridSize + gridSize * .9)
            xOffset += 1
        #text(splitLine[0], 0, (gridSize * yOffset))
        yOffset += 1
        splitLine = readNextLineFromFile()
    delay(20)
    worldTime += 1
    if (splitLine[0] == "new"):
        mapTime += 1
    #noLoop()

    
