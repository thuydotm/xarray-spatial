import numpy as np
import xarray as xr

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
