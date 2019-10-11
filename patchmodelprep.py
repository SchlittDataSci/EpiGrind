import pandas as pd
import geopandas as gp

crsESPG = 4326


### Web GIS Shape Loaders

def getShapeCol(df,index):
    """Gets the shape column for a geo data frame"""
    if index == 'null':
        index = 'geometry'
    return index


def pullWebFile(url):
    """Downloads file and saves it locally, decompressing if needed"""
    filePuller = urllib.URLopener()
    fileName = url.split('/')[-1]
    filePuller.retrieve(url,fileName)
    
    print "Downloading file %s from %s" % (fileName,url.split('://')[-1].split('/')[0])
    zipStyles = ('zip')
    if fileName.split('.')[-1] in zipStyles:
        with zipfile.ZipFile(fileName, "r") as z:
            fileName = fileName.split('.')[0]
            z.extractall()
        print fileName
        listed = listdir(fileName)
        print "Extracted %s to directory %s" % (listed,fileName)
        return ["%s/%s" % (fileName,i) for i in listed]
    
    return [fileName]


def getWebShape(fileRef):
    """Puls a shape file from the internet"""
    if 'http' in fileRef:
        fileRef = pullWebFile(fileRef)
    if type(fileRef) is list:
        fileRef = [i for i in fileRef if i.endswith('.shp')][0]
    return fileRef


def getShapes(fileRef,shapeCol='null'):
    """Extracts shapes universally from a passed fileref"""
    fileRef = getWebShape(fileRef)
    geoDF = gp.GeoDataFrame.from_file(fileRef).to_crs(epsg=crsESPG)
    shapeCol = getShapeCol(geoDF,shapeCol)
    geoSeries = geoDF[shapeCol]
    geoDF[shapeCol] = geoSeries
    return geoDF


### Patch Model Construction Functions

def autoConvert(item):
    try:
        item = float(item)
    except:
        pass
    return item


def divideNumerics(item,scaler):
    try:
        result = item/scaler
        return result
    except:
        return item

    
def checkVal(val,cutoff):
    """Returns true if it is a number above the cutoff or a flag"""
    try:
        val=float(val)
        return val>cutoff
    except:
        return '%(' in val

    
def prepExposures(expRates,patches,cutOff = 0):
    """Generates str exposure script"""
    exposures = dict()

    for toKey in patches:
        exposures[toKey] = {(fromKey,autoConvert(expRates[fromKey][toKey])) for fromKey in patches if checkVal(expRates[fromKey][toKey],cutOff)}
    return str(exposures)


def getEdgeWeights(transFile='Data/Africa_Travel.csv',patches='null',transCol='GravMod'):
    travelMatrix = pd.DataFrame.from_csv(transFile)
    travelPatches = set()
    for pair in travelMatrix['Names'].tolist():
        splitPair = pair.split(' to ')
        travelPatches.update({splitPair[0],splitPair[1]})
    if patches != 'null':
        patchSet = set(patches)
        unionSet = travelPatches.intersection(patchSet)
        diffSet = sorted(list(travelPatches.symmetric_difference(patchSet)))
        print("Prepping sim from union patch set of popfile v travel file.")
        print "Ignoring countries from difference set:\n%s\n\n" % diffSet
        travelPatches = unionSet
        
    maxVal = travelMatrix[transCol].max()
    travelMatrix['WeightedValues'] = [divideNumerics(i,maxVal) for i in travelMatrix[transCol]]
    travelMatrix = travelMatrix.set_index('Names')
    
    return sorted(list(travelPatches)),travelMatrix


def prepCompartments(popFile,countries):
    """Preps str compartments using population data from shapefile"""
    compartments = dict()
    hide = {'compartmentName'}
    toGrab = [item for item in popFile.columns if item not in hide]

    for index, row in popFile.iterrows():
        if row['shapeName'] in countries:
            compartments[index] = {item:row[item] for item in toGrab}
    return str(compartments),compartments.keys()



### Patch Model Generator

def prepScenario(popFile,
                 region=2,
                 selfExposure = 'null',
                 transFile='null',
                 transCol='GravMod',
                 shapeFile = "WorldBorders/TM_WORLD_BORDERS-0.3.shp"):
    """Generates a test scenario from world shapefile and population file"""
    
    shapes = getShapes(shapeFile)
    shapes = shapes[shapes['REGION'] == region]
    countries = shapes['NAME'].tolist()
    
    popData = pd.DataFrame.from_csv(popFile)
    compartmentsStr, patches = prepCompartments(popData,countries)

    if transFile != 'null':
        travelPatches, travelMatrix = getEdgeWeights(transFile,
                                                       patches=patches,
                                                       transCol=transCol)
        shapes = shapes[shapes['NAME'].isin(travelPatches)]
                                         
    shapes = shapes.set_index('NAME')
    centroids = shapes['geometry'].centroid

    expRates = pd.DataFrame(0.,index=travelPatches,columns=travelPatches)
    
    if selfExposure != 'null':
        for i in patches:
            expRates[i][i] = selfExposure
            
    for index,row in travelMatrix.iterrows():
        [i,j] = index.split(' to ')
        value = row['WeightedValues']
        try:
            expRates[i][j] = value
        except:
            pass
    
    exposureStr = prepExposures(expRates,travelPatches)
    
    outStr = 'compartments = %s\nexposures = %s\n#ENABLE CHLOROPLETH\n' % (compartmentsStr,exposureStr)
    outStr = outStr.replace("'%s'" % selfExposure,selfExposure)
    outStr = outStr.replace("'%(",'%(')
    outStr = outStr.replace(")s'",')s')
    return outStr,shapes,countries
