import numpy as np
import xarray as xr
from sortedcontainers import SortedDict

# need to change this value to nan?
NODATA = -1


def equal(a, b):
    return np.isclose(a, b)


class RPolygon:
    # This is a helper class to hold polygons while they are being
    # formed in memory, and to provide services to coalesce a much
    # of edge sections into complete rings.
    def __init__(self):
        self.dfPolyValue = 0
        self.nLastLineUpdated = -1
        self.oMapStrings = SortedDict({})
        self.oMapStartStrings = SortedDict({})
        self.oMapEndStrings = SortedDict({})
        self.iNextStringId = 0
        #     typedef int StringId
        #     typedef std::map< XY, std::pair<StringId, StringId>> MapExtremity
        #
        #     std::map< StringId, std::vector<XY> > oMapStrings{}
        #     MapExtremity oMapStartStrings{}
        #     MapExtremity oMapEndStrings{}

    def findExtremityNot(self, oMap, xy, excludedId):
        if xy == oMap.items()[-1][0]
            return -1
        oIter = oMap.get(xy)
        if oIter[0] != excludedId:
            return oIter[0]
        if oIter[1] != excludedId:
            return oIter[1]
        return -1

    def removeExtremity(self, oMap, xy, id):
        assert xy != oMap.items()[-1][0]

        if oMap[xy][0] == id:
            oMap[xy][0] = oMap[xy][1]
            oMap[xy][1] = -1
            if oMap[xy][0] < 0:
                del oMap[xy]
        elif oMap[xy][1] == id:
            oMap[xy][1] = -1
            assert oMap[xy][0] >= 0
        else:
            assert False
        return

    def insertExtremity(self, oMap, xy, id):
        if (xy != oMap.items()[-1][0]):
            assert oMap[xy][1] == -1
            oMap[xy][0] = id
        else:
            oMap[xy] = (id, -1)
        return

    def merge(self, iBaseString, iSrcString, iDirection):
        anString = self.oMapStrings.get(iSrcString)
        iStart = 1
        iEnd = -1
        if iDirection == 1:
            iEnd = anString.size
        else:
            iStart = anString.size - 2

        self.removeExtremity(self.oMapEndStrings,
                             self.oMapStrings[iBaseString][-1],
                             iBaseString)

        # change size of self.oMapStrings[iBaseString]
        # anBase = self.oMapStrings.get(iBaseString)
        # anBase.reserve(anBase.size() + anString.size() - 1)
        for i in range(iStart, iEnd, iDirection):
            self.oMapStrings[iBaseString].append(anString[i])

        self.removeExtremity(self.oMapStartStrings, anString[0], iSrcString)
        self.removeExtremity(self.oMapEndStrings, anString[-1], iSrcString)
        del self.oMapStrings[iSrcString]
        self.insertExtremity(self.oMapEndStrings,
                             self.oMapStrings[iBaseString][-1],
                             iBaseString)

    def coalesce(self):
        # Iterate over loops starting from the first, trying to merge
        # other segments into them.
        for i in iter(self.oMapStrings):
            thisId = i
            oString = self.oMapStrings.items()[i]
            # Keep trying to merge others strings into our target "base"
            # string while there are matches.
            while True:
                nOtherId = self.findExtremityNot(self.oMapStartStrings,
                                                 oString[-1],
                                                 thisId)
                if nOtherId != -1:
                    self.merge(thisId, nOtherId, 1)
                    continue
                else:
                    nOtherId = self.findExtremityNot(self.oMapEndStrings,
                                                     oString[-1],
                                                     thisId)
                    if nOtherId != -1:
                        self.merge(thisId, nOtherId, -1)
                        continue
                break
            # At this point our loop *should* be closed!
            assert oString[0] == oString[-1]

    def add_segment(self, x1, y1, x2, y2):
        self.nLastLineUpdated = np.max(y1, y2)

        # Is there an existing string ending with this?
        xy1 = (x1, y1)
        xy2 = (x2, y2)

        iExistingString = self.findExtremityNot(self.oMapEndStrings,
                                                xy1, -1)
        if iExistingString >= 0:
            xy1, xy2 = xy2, xy1
        else:
            iExistingString = self.findExtremityNot(self.oMapEndStrings,
                                                    xy2, -1)
        if iExistingString >= 0:
            anString = oMapStrings[iExistingString]
            nSSize = anString.size()

            # We are going to add a segment, but should we just extend
            # an existing segment already going in the right direction?
            nLastLen = np.max(np.abs(anString[nSSize - 2][0] - anString[nSSize - 1][0]),
                              np.abs(anString[nSSize - 2][1] - anString[nSSize - 1][1]))

            self.removeExtremity(self.oMapEndStrings,
                                 anString[-1],
                                 iExistingString)

            if (anString[nSSize - 2][0] - anString[nSSize - 1][0] ==
                (anString[nSSize - 1][0] - xy1[0]) * nLastLen) and \
                    (anString[nSSize - 2][1] - anString[nSSize - 1][1] ==
                     (anString[nSSize - 1][1] - xy1[1]) * nLastLen):
                anString[nSSize - 1] = xy1
            else:
                anString.append(xy1)

            self.insertExtremity(self.oMapEndStrings,
                                 anString[-1],
                                 iExistingString)
            return

        # Create a new string
        self.oMapStrings[self.iNextStringId] = [xy1, xy2]
        self.insertExtremity(self.oMapStartStrings, xy1, self.iNextStringId)
        self.insertExtremity(self.oMapEndStrings, xy2, self.iNextStringId)
        self.iNextStringId += 1


class PolygonEnumerator:
    def __init__(self):
        # TODO: new constructor with param values
        self.nNextPolygonId = 0
        self.nPolyAlloc = 0
        self.connectivity = 0
        self.panPolyIdMap = []
        self.panPolyValue = []

    def new_polygon(self, nValue):
        # Allocate a new polygon id, and reallocate the polygon maps if needed
        nPolyId = self.nNextPolygonId

        # if (nNextPolygonId >= nPolyAlloc)
        #     {
        #         nPolyAlloc = nPolyAlloc * 2 + 20
        #     panPolyIdMap = static_cast < GInt32 * > (
        #         CPLRealloc(panPolyIdMap, nPolyAlloc * sizeof(GInt32)))
        #     panPolyValue = static_cast < DataType * > (
        #         CPLRealloc(panPolyValue, nPolyAlloc * sizeof(DataType)))
        #     }

        self.nNextPolygonId += 1
        self.panPolyIdMap[nPolyId] = nPolyId
        self.panPolyValue[nPolyId] = nValue

    def merge_polygon(self, nSrcId, nDstIdInit):
        # Update the polygon map to indicate the merger of two polygons.
        # Figure out the final dest id.
        nDstIdFinal = nDstIdInit
        while (self.panPolyIdMap[nDstIdFinal] != nDstIdFinal):
            nDstIdFinal = self.panPolyIdMap[nDstIdFinal]

        # Map the whole intermediate chain to it.
        nDstIdCur = nDstIdInit
        while (self.panPolyIdMap[nDstIdCur] != nDstIdCur):
            nNextDstId = self.panPolyIdMap[nDstIdCur]
            self.panPolyIdMap[nDstIdCur] = nDstIdFinal
            nDstIdCur = nNextDstId

        # And map the whole source chain to it too (can be done in one pass).
        while (self.panPolyIdMap[nSrcId] != nSrcId):
            nNextSrcId = self.panPolyIdMap[nSrcId]
            self.panPolyIdMap[nSrcId] = nDstIdFinal
            nSrcId = nNextSrcId

        self.panPolyIdMap[nSrcId] = nDstIdFinal
        return

    def process_line(self, panLastLineVal, panThisLineVal,
                     panLastLineId, panThisLineId, nXSize):
        # special case for the first line
        if panLastLineVal is None:
            for i in range(nXSize):
                if panThisLineVal[i] == NODATA:
                    # check if nodata cell
                    panThisLineId[i] = -1
                elif i == 0 or not equal(panThisLineVal[i],
                                         panThisLineVal[i - 1]):
                    panThisLineId[i] = self.new_polygon(panThisLineVal[i])
                else:
                    panThisLineId[i] = panThisLineId[i - 1]
        # Process each pixel comparing to the previous pixel, and
        # to the last line.
        for i in range(nXSize):
            if panThisLineVal[i] == NODATA:
                # check if nodata cell
                panThisLineId[i] = -1

            elif i > 0 and equal(panThisLineVal[i], panThisLineVal[i - 1]):
                panThisLineId[i] = panThisLineId[i - 1]

                if equal(panLastLineVal[i], panThisLineVal[i]) and \
                        self.panPolyIdMap[panLastLineId[i]] != \
                        self.panPolyIdMap[panThisLineId[i]]:
                    self.merge_polygon(panLastLineId[i], panThisLineId[i])

                if self.connectivity == 8 and equal(panLastLineVal[i - 1],
                                                    panThisLineVal[i]) \
                        and self.panPolyIdMap[panLastLineId[i - 1]] != \
                        self.panPolyIdMap[panThisLineId[i]]:
                    self.merge_polygon(panLastLineId[i - 1], panThisLineId[i])

                if self.connectivity == 8 and i < nXSize - 1 and \
                        equal(panLastLineVal[i + 1], panThisLineVal[i]) and \
                        self.panPolyIdMap[panLastLineId[i + 1]] != \
                        self.panPolyIdMap[panThisLineId[i]]:
                    self.merge_polygon(panLastLineId[i + 1], panThisLineId[i])

            elif equal(panLastLineVal[i], panThisLineVal[i]):
                panThisLineId[i] = panLastLineId[i]

            elif i > 0 and self.connectivity == 8 and \
                    equal(panLastLineVal[i - 1], panThisLineVal[i]):
                panThisLineId[i] = panLastLineId[i - 1]
                if i < nXSize - 1 and \
                        equal(panLastLineVal[i + 1], panThisLineVal[i]) and \
                        self.panPolyIdMap[panLastLineId[i + 1]] != \
                        self.panPolyIdMap[panThisLineId[i]]:
                    self.merge_polygon(panLastLineId[i + 1], panThisLineId[i])

            elif i < nXSize - 1 and self.connectivity == 8 and \
                    equal(panLastLineVal[i + 1], panThisLineVal[i]):
                panThisLineId[i] = panLastLineId[i + 1]

            else:
                panThisLineId[i] = self.new_polygon(panThisLineVal[i])
        return

    def complete_merges(self):
        # Make a pass through the maps, ensuring every polygon id
        # points to the final id it should use, not an intermediate value.
        nFinalPolyCount = 0
        for iPoly in range(self.nNextPolygonId):
            # Figure out the final id.
            nId = self.panPolyIdMap[iPoly]
            while (nId != self.panPolyIdMap[nId]):
                nId = self.panPolyIdMap[nId]
                
            # Then map the whole intermediate chain to it.
            nIdCur = self.panPolyIdMap[iPoly]
            self.panPolyIdMap[iPoly] = nId
            while (nIdCur != self.panPolyIdMap[nIdCur]):
                nNextId = self.panPolyIdMap[nIdCur]
                self.panPolyIdMap[nIdCur] = nId
                nIdCur = nNextId

            if self.panPolyIdMap[iPoly] == iPoly:
                nFinalPolyCount += 1
                
        return


def swap(arr1, arr2):
    arr1, arr2 = arr2, arr1
    return 


def add_edges(panThisLineId, panLastLineId,
              panPolyIdMap, panPolyValue,
              papoPoly, iX, iY):

    # Examine one pixel and compare to its neighbour above
    # (previous) and right.  If they are different polygon ids
    # then add the pixel edge to this polygon and the one on the
    # other side of the edge.   

    # TODO(schwehr): Simplify these three vars.
    nThisId = panThisLineId[iX]
    if nThisId != -1:
        nThisId = panPolyIdMap[nThisId]
    nRightId = panThisLineId[iX+1]
    if nRightId != -1:
        nRightId = panPolyIdMap[nRightId]
    nPreviousId = panLastLineId[iX]
    if nPreviousId != -1:
        nPreviousId = panPolyIdMap[nPreviousId]

    iXReal = iX - 1

    if nThisId != nPreviousId:
        if nThisId != -1:
            if papoPoly[nThisId] == nullptr:
                papoPoly[nThisId] = RPolygon(panPolyValue[nThisId])
            papoPoly[nThisId].add_segment(iXReal, iY, iXReal+1, iY)

        if nPreviousId != -1:
            if papoPoly[nPreviousId] == nullptr:
                papoPoly[nPreviousId] = RPolygon(panPolyValue[nPreviousId])
            papoPoly[nPreviousId].add_segment(iXReal, iY, iXReal+1, iY)

    if nThisId != nRightId:
        if nThisId != -1:
            if papoPoly[nThisId] == nullptr:
                papoPoly[nThisId] = RPolygon(panPolyValue[nThisId])
            papoPoly[nThisId].add_segment(iXReal+1, iY, iXReal+1, iY+1)

        if nRightId != -1:
            if papoPoly[nRightId] == nullptr:
                papoPoly[nRightId] = RPolygon(panPolyValue[nRightId])
            papoPoly[nRightId].add_segment(iXReal+1, iY, iXReal+1, iY+1)
    return


def emit_polygon_to_layer(hOutLayer, iPixValField,
                          poRPoly, padfGeoTransform):
    # Turn bits of lines into coherent rings.
    poRPoly.coalesce()

    # Create the polygon geometry.
    OGRGeometryH hPolygon = OGR_G_CreateGeometry( wkbPolygon );




def polygonize(raster, mask, connectivity):
    """

    """
    height, width = raster.shape
    panLastLineVal = np.zeros(width + 2)
    panThisLineVal = np.zeros(width + 2)
    panLastLineId = np.zeros(width + 2)
    panThisLineId = np.zeros(width + 2)

    paby_mask_line = None
    if mask is not None:
        paby_mask_line = np.zeros(width)

    # Get the geotransform, if there is one, so we can convert the
    #   vectors into georeferenced coordinates.
    if 'affine' in raster.dims:
        affine_coefficients = raster['affine']
    else:
        affine_coefficients = [0, 1, 0, 0, 0, 1]

    # The first pass over the raster is only used to build up the
    # polygon id map so we will know in advance what polygons are
    # what on the second pass.
    oFirstEnum = PolygonEnumerator()

    for iY in range(height):
        # read line by line
        panThisLineVal = raster.data[iY]
        if iY == 0:
            oFirstEnum.process_line(None, panThisLineVal,
                                    None, panThisLineId, width)
        else:
            oFirstEnum.process_line(panLastLineVal, panThisLineVal,
                                    panLastLineId, panThisLineId,
                                    width)
        # Swap lines
        swap(panLastLineVal, panThisLineVal)
        swap(panLastLineId, panThisLineId)
        
    # Make a pass through the maps, ensuring every polygon id
    # points to the final id it should use, not an intermediate value.
    oFirstEnum.complete_merges()
    
    # Initialize ids to -1 to serve as a nodata value for the
    # previous line, and past the beginning and end of the scanlines.
    nXSize = width
    panThisLineId[0] = -1
    panThisLineId[nXSize + 1] = -1
    for iX in range(nXSize+2):
        panLastLineId[iX] = -1

    # We will use a new enumerator for the second pass primarily
    # so we can preserve the first pass map.
    oSecondEnum = PolygonEnumerator()
    # RPolygon ** papoPoly = static_cast < RPolygon ** > (
    #     CPLCalloc(sizeof(RPolygon *), oFirstEnum.nNextPolygonId))
    papoPoly = [RPolygon() for i in range(oFirstEnum.nNextPolygonId)]
    # Second pass during which we will actually collect polygon
    # edges as geometries

    for iY in range(nYSize+1):
        # Read the image data.                                            */
        if iY < nYSize:
            panThisLineVal = hSrcBand[iY]
            # if eErr == CE_None && hMaskBand != nullptr
            #     eErr = GPMaskImageData(hMaskBand, pabyMaskLine, iY, nXSize,
            #                             panThisLineVal
        # Determine what polygon the various pixels belong to (redoing
        # the same thing done in the first pass above).
        if iY == nYSize:
            for iX in range(nXSize+2):
                panThisLineId[iX] = -1
        elif iY == 0:
            oSecondEnum.process_line(nullptr, panThisLineVal, nullptr,
                                     panThisLineId+1, nXSize)
        else:
            oSecondEnum.process_line(panLastLineVal, panThisLineVal,
                                     panLastLineId+1,  panThisLineId+1,
                                     nXSize)

        # Add polygon edges to our polygon list for the pixel
        # boundaries within and above this line.
        for iX in range(nXSize+1):
            AddEdges(panThisLineId, panLastLineId,
                     oFirstEnum.panPolyIdMap, oFirstEnum.panPolyValue,
                     papoPoly, iX, iY)

        # Periodically we scan out polygons and write out those that
        # haven't been added to on the last line as we can be sure
        # they are complete.
        if iY % 8 == 7:
            for iX in range(oSecondEnum.nNextPolygonId):
                if papoPoly[iX] is not None and papoPoly[iX].nLastLineUpdated < iY-1:
                    eErr =
                        EmitPolygonToLayer(hOutLayer, iPixValField,
                                            papoPoly[iX], adfGeoTransform

                    delete papoPoly[iX]
                    papoPoly[iX] = nullptr

        # Swap pixel value, and polygon id lines to be ready for the      */
        # next line.                                                      */
        swap(panLastLineVal, panThisLineVal)
        swap(panLastLineId, panThisLineId)

    # Make a cleanup pass for all unflushed polygons.                 */
    for iX in range(oSecondEnum.nNextPolygonId):
        if papoPoly[iX] is not None:
            eErr = EmitPolygonToLayer(hOutLayer, iPixValField,
                                       papoPoly[iX], adfGeoTransform

            del papoPoly[iX]
            papoPoly[iX] = nullptr
            
    return