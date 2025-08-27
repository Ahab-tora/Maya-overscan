#--- --- Imports --- ---#
import maya.OpenMaya as om
import maya.cmds as cmds
from math import sqrt,atan,pi

#--- --- --- ---#

'''
Example use:
'''


class OverscanCamera():


    def __init__(self,newCamName:str,
                 firstFrame:int, lastFrame:int,
                 horizontalResolution:int, verticalResolution:int,
                 inputCam:str):

        self.newCamName:str = newCamName
        self.frameRange = [firstFrame,lastFrame]
        self.horizontalResolution:int = horizontalResolution
        self.verticalResolution:int = verticalResolution
        self.inputCam:str = inputCam
        self.farClip:int = -100000 #if there is no collision , the collision point will be at farClip on all axis


    def getCamData(self,cam:str) -> dict:

        '''
        gets data from a cam: horizontal and vertical aperture, and focal length
        '''

        horizontalAperture = cmds.getAttr(cam + '.horizontalFilmAperture')*2.54*10 #inches to cm ? or meters? need to check
        verticalAperture = cmds.getAttr(cam + '.verticalFilmAperture')*2.54*10
        focalLength  =  cmds.getAttr(cam + '.focalLength')

        dataDict = {'horizontalAperture': horizontalAperture,
                    'verticalAperture' : verticalAperture,
                    'focalLength' : focalLength}

        return dataDict


    def createProjectionPoints(self,
                            cam:str,
                            rows:int=3,columns:int=3) -> list:

        '''
        creates a grid of locators parented to the camera used for the raycasting
        returns a list of the locs
        '''

        camDict = self.getCamData(cam = cam)

        rows = rows - 1
        columns = columns - 1

        viewPlane = cmds.polyPlane(n = 'cam_viewPlane',sx = rows,sy = columns )[0]
        cmds.parent(viewPlane, cam, r=True)
        
        cmds.setAttr(viewPlane + '.rotateX',90)
        cmds.setAttr(viewPlane + '.scaleX', camDict['horizontalAperture'] ) #creates a viewplane "similar to a lens" to scatter the locs on
        cmds.setAttr(viewPlane + '.scaleZ', camDict['verticalAperture'])
        cmds.setAttr(viewPlane + '.translateZ', -(camDict['focalLength']))

        projectionPoints = []

        for i,vtx in enumerate(cmds.ls(viewPlane + '.vtx[*]', flatten = True)):

            translation = cmds.xform(vtx,q=True,t=True,ws=True)
            locName = 'projectionPoint' + str(i)
            loc = cmds.spaceLocator(n=locName ,p=(0,0,0))[0]
            projectionPoints.append(loc)
            
            for axis,val in zip(['X','Y','Z'],translation):
                cmds.setAttr(loc + '.translate' + axis,val)
                
            cmds.parent(loc,cam)

        cmds.delete(viewPlane) #not needed anymore

        return projectionPoints


    def getMeshes(self):
        'returns all the meshes shapes with the visibility ON'

        meshesShapes = []
        allMeshesShapes = cmds.ls(type='mesh')

        for meshShape in allMeshesShapes:
            transformNode = cmds.listRelatives(meshShape,p=True)[0]
            visState = cmds.getAttr(transformNode + '.visibility')
            if visState:
                meshesShapes.append(meshShape)

        return meshesShapes


    def createAvgCam(self):

        '''
        returns a group 'avgCamSpace'  at the average position of the inputCam
        '''

        def buildCamPosDict(cam:str,frameRange:list[int]) -> dict:
            'build a dict with the translation/rotation (worldSpace) in a framerange'

            furthestCam = 0 #check that
            camPosition = {}

            for frame in range(frameRange[0],frameRange[1]+1):
                cmds.currentTime(frame)

                camTX,camTY,camTZ  = cmds.xform(cam,q=True,t=True,ws=True)
                camRX,camRY,camRZ = cmds.xform(cam,q=True,ro=True,ws=True)

                if furthestCam < camTZ:
                    furthestCam = camTZ

                camPosition[frame] = {}
                camPosition[frame]['translation'] = {'X':camTX,'Y':camTY,'Z':camTZ}
                camPosition[frame]['rotation'] = {'X':camRX,'Y':camRY,'Z':camRZ}

            return camPosition,furthestCam
        
        def computeAvgCam(camPosition:dict):
            'takes an a dict of the cam position through multiples frames and returns the average position of the cam'

            valsTX = []
            valsTY = []
            valsTZ = []

            valsRX = []
            valsRY = []
            valsRZ = []

            for frame in camPosition.keys():
                valsTX.append(camPosition[frame]['translation']['X'])
                valsTY.append(camPosition[frame]['translation']['Y'])
                valsTZ.append(camPosition[frame]['translation']['Z'])

                valsRX.append(camPosition[frame]['rotation']['X'])
                valsRY.append(camPosition[frame]['rotation']['Y'])
                valsRZ.append(camPosition[frame]['rotation']['Z'])

            avgTX = sum(valsTX)/len(valsTX)
            avgTY = sum(valsTY)/len(valsTY)
            avgTZ = sum(valsTZ)/len(valsTZ)

            avgRX = sum(valsRX)/len(valsRX)
            avgRY = sum(valsRY)/len(valsRY)
            avgRZ = sum(valsRZ)/len(valsRZ)

            avgCamTranslation = {'X':avgTX,'Y':avgTY,'Z':avgTZ}
            avgCamRotation = {'X':avgRX,'Y':avgRY,'Z':avgRZ}

            return avgCamTranslation,avgCamRotation
        
        camPosition,furthestCam = buildCamPosDict(cam=self.inputCam, frameRange=self.frameRange)
        avgCamTranslation,avgCamRotation = computeAvgCam(camPosition=camPosition)
        avgCamSpace = cmds.group(n='avgCamSpace',em=True)

        for axis in ['X','Y','Z']:

            cmds.setAttr(avgCamSpace + '.translate' + axis,avgCamTranslation[axis])
            cmds.setAttr(avgCamSpace + '.rotate' + axis,avgCamRotation[axis])

        return avgCamSpace


    def closestIntersection(self,
                             source:tuple,vector:tuple,
                             meshesShapes:list,name:str):
        

        '''
        raycast between source and vector through the meshes shapes, returns a locator at the closest intersection
        '''

        nearestIntersectionDistance = 999999999
        intersectionPoint = (None,None,None)

        for meshShape in meshesShapes:

            try:
                data = self.ray(mesh=meshShape, source=source,
                           vector = vector, farClip = self.farClip)
                
            except Exception as e:
                print(f'Error while raycasting: {e}')

            intersectionSourceDistance = self.distanceBetweenTwoPoints3D(point1=(data['X'],data['Y'],data['Z'],),
                                                point2=source)

            if intersectionSourceDistance <= nearestIntersectionDistance:
                nearestIntersectionDistance = intersectionSourceDistance
                intersectionPoint = (data['X'],data['Y'],data['Z'],)

        projectedLoc = cmds.spaceLocator(n = name, p=(0,0,0))[0]
        cmds.xform(projectedLoc, ws=True, t=(intersectionPoint[0], intersectionPoint[1], intersectionPoint[2]))

        return projectedLoc


    def angleBetweenPoints(self,sourcePos:tuple,pointPos:tuple,sqrPoint:int):

        '''
        using trigo, check what the angle is between two points, sqrPoints should be aligned on two axis with the source, and on one axis with pointPos
        returns an angle (degrees)
        '''
        adjDistance = self.distanceBetweenTwoPoints2D(point1=pointPos,point2=(sourcePos[0],sqrPoint)) #adjacent to the angle
        opDistance = self.distanceBetweenTwoPoints2D(point1=sourcePos,point2=pointPos) #opposite to the angle
            
        radiantVal = atan(adjDistance/opDistance) #returns the value as radiant
        degreesVal = radiantVal * 180 / pi #translating radiant to degrees
        if pointPos[1]:
            degreesVal = degreesVal *(-1)
        return degreesVal


    def midPoints(self,points:list[tuple])->tuple:
        '''
        Takes a list of 2 points in a 3D space (X,Y,Z) and returns the point in the middle of them
        '''
        pointsNb = len(points)
        xVals = []
        yVals = []
        zVals = []
        for point in points:
            xVals.append(point[0])
            yVals.append(point[1])
            zVals.append(point[2])

        averageX = sum(xVals)/pointsNb
        averageY = sum(yVals)/pointsNb
        averageZ = sum(zVals)/pointsNb

        return(averageX,averageY,averageZ,)


    def distanceBetweenTwoPoints3D(self,point1:tuple,point2:tuple):
        '''returns the distance between two points'''
        dX = point1[0] - point2[0]
        dY = point1[1] - point2[1]
        dZ = point1[2] - point2[2]
        return sqrt(dX * dX + dY * dY + dZ * dZ)


    def distanceBetweenTwoPoints2D(self,point1:tuple,point2:tuple) -> tuple :
        'returns the distance between two points in a 2D space (x,y)'
        x = point1[0] - point2[0]
        y = point1[1] - point2[1]
        return(x**2 + y**2)


    def ray(self,mesh,source:tuple,vector:tuple,farClip:int,printMesh:bool=False) -> dict:
        'casts a ray from a source through a vector to a mesh to check for any intersection'

        cmds.select(cl=True)
        item = om.MDagPath()

        selList = om.MSelectionList()
        selList.add(mesh)
        
        selList.getDagPath(0, item)
        
        item.extendToShape()
        fnMesh = om.MFnMesh(item)

        raySource = om.MFloatPoint(source[0], source[1], source[2], 1.0)

        rayDirection = om.MFloatVector(vector[0]-source[0], vector[1]-source[1], vector[2]-source[2])

        faceIds = None
        triIds = None
        idsSorted = True

        space = om.MSpace.kWorld

        maxParam = 10001
        testBothDirections = False
        accelParams = None
        sortHits = True

        hitPoints = om.MFloatPointArray()

        hitRayParams =  om.MFloatArray()
        hitFaces = om.MIntArray()
        hitTriangles = None
        hitBary1s = None
        hitBary2s = None
        tolerance = 0.0001

        
        hit = fnMesh.allIntersections(raySource, rayDirection, faceIds, triIds, idsSorted,
                                    space, maxParam, testBothDirections,
                                    accelParams, sortHits, hitPoints,
                                    hitRayParams, hitFaces, hitTriangles, hitBary1s,
                                    hitBary2s, tolerance)
        
        if hit:
            dataDict = {
                        'X':hitPoints[0][0],
                        'Y':hitPoints[0][1],
                        'Z':hitPoints[0].z,
                        'hit':True,
                        'mesh':mesh
                        }
            if printMesh:
                print(f"Collision on : {mesh}")

            
        else:
            dataDict = {
                        'X':farClip,
                        'Y':farClip,
                        'Z':farClip,
                        'hit':False,
                        'mesh':None
                        }
            
        return dataDict


    def overscan(self):
        '''
        main function
        '''

        meshesShapes = self.getMeshes() 
        
        projectionPoints = self.createProjectionPoints(cam = self.inputCam) 

        avgCamSpace = self.createAvgCam()

        farthestYAngle = {'angle':0, 'point':None}
        farthestZAngle = {'angle':0, 'point':None}
        closestZ = self.farClip

        projectedLocs = []
        for frame in range(self.frameRange[0],self.frameRange[1]+1): #for each frame, create a loc if there is an intersection between a mesh and

            cmds.currentTime(frame)
            source = cmds.xform(self.inputCam ,q=True, t=True, ws=True)

            for projectionPoint in projectionPoints:

                vector = cmds.xform(projectionPoint, q=True, t=True, ws=True)


                projectecLoc = self.closestIntersection(source = source,vector = vector,
                                                        meshesShapes= meshesShapes,name = projectionPoint + '_' + str(frame))
                projectedLocs.append(projectecLoc)
                cmds.parent(projectecLoc,avgCamSpace)
                
                origin = (0,0,0)
                projectecLocPos = cmds.getAttr(projectecLoc + '.translate')[0]

                angleY = self.angleBetweenPoints(sourcePos=(origin[1],origin[2]),
                    pointPos=(projectecLocPos[1],projectecLocPos[2]),sqrPoint= -abs(projectecLocPos[2])) #checking if the angle between the cam and the locator is the largest on a Y rotation 
                
                angleZ = self.angleBetweenPoints(sourcePos=(origin[0],origin[2]),
                    pointPos=(projectecLocPos[0],projectecLocPos[2]),sqrPoint=projectecLocPos[2]) #same on Z rotation


                if abs(farthestYAngle['angle']) < abs(angleY):
                    farthestYAngle['angle'] = angleY         
                    farthestYAngle['point'] = projectecLoc
                    if projectecLocPos[2] > closestZ:
                        closestZ = projectecLocPos[2]

                if abs(farthestZAngle['angle']) < abs(angleZ):
                    farthestZAngle['angle'] = angleZ
                    farthestZAngle['point'] = projectecLoc
                    if projectecLocPos[2] > closestZ:
                        closestZ = projectecLocPos[2]

        farthestXpos = cmds.getAttr(farthestZAngle['point'] + '.translateX')
        farthestYAnglepos = cmds.getAttr(farthestYAngle['point'] + '.translateY')

        #4 points used to get the new aperture : TL,TR,BL,BR
        TL = cmds.spaceLocator(n='TL', p=(0,0,0))[0]
        TLpos = (-abs(farthestXpos),abs(farthestYAnglepos),closestZ)
        cmds.parent(TL,avgCamSpace,r=True)
        cmds.xform(TL,t=TLpos)

        TR = cmds.spaceLocator(n='TR', p=(0,0,0))[0]
        TRpos = (abs(farthestXpos),abs(farthestYAnglepos),closestZ)
        cmds.parent(TR,avgCamSpace)
        cmds.xform(TR,t=TRpos)

        BL = cmds.spaceLocator(n='BL', p=(0,0,0))[0]
        BLpos = (-abs(farthestXpos),-abs(farthestYAnglepos),closestZ)
        cmds.parent(BL,avgCamSpace)
        cmds.xform(BL,t=BLpos)

        BR = cmds.spaceLocator(n='BR', p=(0,0,0))[0]
        BRpos = (abs(farthestXpos),-abs(farthestYAnglepos),closestZ)
        cmds.parent(BR,avgCamSpace)
        cmds.xform(BR,t=BRpos)

        outputCamName = self.newCamName
        outputCam = cmds.camera(n=outputCamName)[0]
        cmds.parent(outputCam,avgCamSpace,r=True)




        middlePoint = self.midPoints(points=[TLpos,BRpos])

        distCamMiddlePoint = self.distanceBetweenTwoPoints3D(point1=(0,0,0),
                                point2=middlePoint)


        horizontalDistance = self.distanceBetweenTwoPoints3D(point1=TLpos,
                                                point2=TRpos)


        verticalDistance = self.distanceBetweenTwoPoints3D(point1=TRpos,
                                                point2=BRpos)



        outputCamFocalLength = cmds.getAttr(outputCam + '.focalLength')
        focalFactor = outputCamFocalLength/distCamMiddlePoint

        outputCamVerticalAperture = (verticalDistance * focalFactor)/10/2.54 #dividing by 2.54 to convert to inches, dividing by 10 is needed but no idea why
        outputCamHorizontalAperture = (horizontalDistance * focalFactor)/10/2.54 

        cmds.setAttr(outputCam + '.verticalFilmAperture',outputCamVerticalAperture) 
        cmds.setAttr(outputCam + '.horizontalFilmAperture',outputCamHorizontalAperture)
        cmds.setAttr(outputCam+'.displayCameraFrustum',1)

        inputCamVerticalAperture = cmds.getAttr(self.inputCam + '.verticalFilmAperture')
        inputCamHorizontalAperture = cmds.getAttr(self.inputCam + '.horizontalFilmAperture')
        inputCamFocalLength = cmds.getAttr(self.inputCam + '.focalLength')

        inputCamVerticalApertureNormalized = inputCamVerticalAperture / inputCamFocalLength
        inputCamHorizontalApertureNormalized = inputCamHorizontalAperture / inputCamFocalLength

        outputCamVerticalApertureNormalized = outputCamVerticalAperture / outputCamFocalLength
        outputCamHorizontalApertureNormalized = outputCamHorizontalAperture / outputCamFocalLength


        newVerticalResolution = round((outputCamVerticalApertureNormalized / inputCamVerticalApertureNormalized) * self.verticalResolution)
        newHorizontalResolution = round((outputCamHorizontalApertureNormalized / inputCamHorizontalApertureNormalized ) * self.horizontalResolution)

        resolutionGRP = cmds.group(n=f'_{newHorizontalResolution} by {newVerticalResolution}_',em=True)
        cmds.parent(resolutionGRP,avgCamSpace)
        print(f'The new resolution is {newHorizontalResolution}*{newVerticalResolution}')

        #deleting what is not needed anymore
        for projectionPoint in projectionPoints:
            cmds.delete(projectionPoint)

        for projectedLoc in projectedLocs:
            cmds.delete(projectedLoc)

        for cornerLoc in TL,TR,BL,BR:
            cmds.delete(cornerLoc)



