function computeTriangulation(points) {
	ini(points);
	triangulate();
	outputTriangles = getOutputTriangles();
	return outputTriangles;
}

function triangulate(){
	for(var i=3; i < vertexs.length; i+=1){
		if(i != fixedPoint){
			addPTT(i)
		}
	}
}

function createEdge(index,vS,vE,fL,fR,eB,eN){
	edges[index] = {
	'vertexStart': vS,
	'vertexEnd': vE,
	'faceLeft': fL,
	'faceRight': fR,
	'edgeBefore': eB,
	'edgeNext': eN,
	};
}

function ini(points){
	nPoints = points.length;
	faces = [nPoints - 2];
	edges = new Array(3 * nPoints - 6);
	vertexs = new Array(nPoints);
	for(var i = 0; i < points.length; i++){
		vertexs[i] = {Coordinates: points[i]}
	}
	countF = 0;
	countE = 0;
	createEdge(0,0,1,0,-1,2,1);
	createEdge(1,1,2,0,-1,0,2);
	createEdge(2,2,0,0,-1,1,0);
	countE += 3;
	faces[0] = 0;
	countF += 1;
	vertexs[0].eI = 0;
	vertexs[1].eI = 1;
	vertexs[2].eI = 2;
	fixedPoint = chooseFixedPointIndex(points);
	fixedCoordinates = vertexs[fixedPoint].Coordinates;
	addPointToFace(fixedPoint, 0);
	countT = 4;
}

function getSVFERL(eI, fI){//
	if (edges[eI].faceLeft == fI){
		return edges[eI].vertexStart;
	}
	else if(edges[eI].faceRight == fI){
		return edges[eI].vertexEnd;
	}
}

function newOne(eI, pfI, nfI){//
	if (edges[eI].faceLeft == pfI){
		edges[eI].faceLeft = nfI
	}
	else {
		edges[eI].faceRight = nfI
	}
}

function replaceEdge(eI, pI, nI){//
	if (edges[eI].edgeBefore == pI){
		edges[eI].edgeBefore = nI
	}
	else {
		edges[eI].edgeNext = nI
	}
}

function getVertexCoordinates(vI){
	return vertexs[vI].Coordinates
}

function getOtherFaceFromEdge(eI, firstFace){//
	var rightFace = edges[eI].faceRight;
	if (rightFace === firstFace){
		return edges[eI].faceLeft;
	}
	return rightFace;
}

function getEEF(fI){//
	var fEF = faces[fI];
	var sEF, tEF;
	if (edges[fEF].faceLeft == fI){
		sEF = edges[fEF].edgeBefore;
		if (edges[fEF].vertexStart == edges[sEF].vertexEnd){
			tEF = edges[sEF].edgeBefore;
		}
		else{
			tEF = edges[sEF].edgeNext;
		}
	}

	else { 
		sEF = edges[fEF].edgeNext;
		if (edges[sEF].vertexStart == edges[fEF].vertexEnd){
			tEF = edges[sEF].edgeNext;
		}
		else{
			tEF = edges[sEF].edgeBefore;
		}
	}
	return [tEF, sEF, fEF];
}

function getVIF(fI) {/////////$$$$$$$$·"·$·"$"·$"No esta bien
	var vertices = [];
	var fEF = faces[fI];
	vertices.push(edges[fEF].vertexStart, edges[fEF].vertexEnd);
	var sEF;

	if (fEF.faceLeft === fI) {
		sEF = fEF.edgeBefore;
	}
	else {
		sEF = edges[fEF].edgeNext;
	}
	vertices.push(edges[sEF].vertexStart, edges[sEF].vertexEnd);
	return Array.from(new Set(vertices))
}

function getFAV(vI){//
	var facesAroundVertex = [];
	var firstEdge = vertexs[vI].eI;
	var edge = firstEdge;
	var tmpFace = null;
	do{
		if (edges[edge].vertexStart == vI){
			tmpFace = edges[edge].faceLeft;
			edge = edges[edge].edgeBefore;
		}
		else {
			tmpFace = edges[edge].faceRight;
			edge = edges[edge].edgeNext;
		}
		facesAroundVertex.push(tmpFace);
	} while (edge !== firstEdge);
	return facesAroundVertex
}

function getOutputTriangles() {//
	var outputTriangles = [];
	for(var i = 0; i<faces.length; i++){
		var faceEdge = faces[i];
		if (faceEdge !== undefined){
			outputTriangles.push(getVIF(i))
		}
	}
	return outputTriangles
}

function chooseFixedPointIndex(points){//
	var minimumXCoordinate = points[0].x;
	var maximumXCoordinate = points[0].x;
	var minimumYCoordinate = points[0].y;
	var maximumYCoordinate = points[0].y;
	for(var i=3; i < points.length; i++){ 
		minimumXCoordinate = Math.min(minimumXCoordinate, points[i].x);
		maximumXCoordinate = Math.max(maximumXCoordinate, points[i].x);
		minimumYCoordinate = Math.min(minimumYCoordinate, points[i].y);
		maximumYCoordinate = Math.max(maximumYCoordinate, points[i].y);
	}
	var xSpan = (maximumXCoordinate - minimumXCoordinate) / 5;
	var ySpan = (maximumYCoordinate - minimumYCoordinate) / 5;
	var xValuesInTheMiddleStart = (minimumXCoordinate + maximumXCoordinate) / 2 - xSpan;
	var xValuesInTheMiddleEnd = (minimumXCoordinate + maximumXCoordinate) / 2 + xSpan;
	var yValuesInTheMiddleStart = (minimumYCoordinate + maximumYCoordinate) / 2 - ySpan;
	var yValuesInTheMiddleEnd = (minimumYCoordinate + maximumYCoordinate) / 2 + ySpan;
	var matchingX;
	for(var i=3; i < points.length; i++){
		if ( (xValuesInTheMiddleStart < points[i].x) && (points[i].x  < xValuesInTheMiddleEnd)){
			matchingX = i;
			if ( (yValuesInTheMiddleStart < points[i].y) && (points[i].y  < yValuesInTheMiddleEnd) ) {
				return i;
			}
		}

	}
	if (matchingX)
		return matchingX;
	return 3;
}

function getEDP(edgesIndices){//
	for(var i=0; i < edgesIndices.length; i++){
		var eI = edgesIndices[i];
		if (edges[eI].vertexStart !== fixedPoint && edges[eI].vertexEnd !== fixedPoint){
			return eI
		}
	}
}


function addPointOnBoundary(eS, nPI){//
	var startVertex = edges[eS].vertexStart;
	var endVertex = edges[eS].vertexEnd;
	var nextEdge = edges[eS].edgeNext;
	var beforeEdge = edges[eS].edgeBefore;
	var leftFace = edges[eS].faceLeft;
	var rightFace = edges[eS].faceRight;
	var tELF = getEEF(leftFace).filter(
		edge => edge !== eS && edge !== beforeEdge
	)[0];
	var thirdEdgeInRightFace = getEEF(rightFace).filter(
		edge => edge !== eS && edge !== nextEdge
	)[0];
	var vALF = getVIF(leftFace);
	var vARF = getVIF(rightFace);
	var thirdVertexFromLeftFace = vALF.filter(
		vertex => vertex !== startVertex && vertex !== endVertex
	)[0];
	var thirdVertexFromRightFace = vARF.filter(
		vertex => vertex !== startVertex && vertex !== endVertex
	)[0];
	var countEBeforeAdding = countE;
	var countFBeforeAdding = countF;
	countE += 3;
	countF += 2;
	vertexs[nPI].eI = countEBeforeAdding;
	faces[countFBeforeAdding] = countEBeforeAdding + 2;
	faces[countFBeforeAdding + 1] = countEBeforeAdding + 2;
	createEdge(countEBeforeAdding,nPI,thirdVertexFromLeftFace,countFBeforeAdding,leftFace,tELF,countEBeforeAdding);
	createEdge(countEBeforeAdding + 1,nPI,thirdVertexFromRightFace,leftFace,countFBeforeAdding + 1,thirdEdgeInRightFace,eS);
	createEdge(countEBeforeAdding + 2,startVertex,nPI,countFBeforeAdding,countFBeforeAdding + 1,countEBeforeAdding + 1,beforeEdge);
	edges[eS].vertexStart = nPI;
	edges[eS].edgeBefore = countEBeforeAdding;
	newOne(beforeEdge, leftFace, countFBeforeAdding);
	newOne(thirdEdgeInRightFace, rightFace, countFBeforeAdding + 1);
}

function addPTT(nPI){
	var newCoordinates = vertexs[nPI].Coordinates;
	var facesAroundFixedVertex = getFAV(fixedPoint);
	for(var i = 0; i < facesAroundFixedVertex.length; i++){
		var faceAroundFixedVertex = facesAroundFixedVertex[i];
		var verticesIndicesAroundFace = getVIF(faceAroundFixedVertex);
		var verticesCoordinatesAroundFace = [];
		for(var j = 0; j < verticesIndicesAroundFace.length; ++j){
			verticesCoordinatesAroundFace[j] = vertexs[verticesIndicesAroundFace[j]].Coordinates;
		}
		if (isInsideTriangle(verticesCoordinatesAroundFace, newCoordinates)){
			addPointToFace(nPI, faceAroundFixedVertex);
			countT++;
		}
		else {
			var edgesAroundFace = getEEF(faceAroundFixedVertex);
			var edgeNotPointingToFixedVertex = getEDP(edgesAroundFace);
			var startVertexCoordinates = vertexs[edges[edgeNotPointingToFixedVertex].vertexStart].Coordinates;;
			var endVertexCoordinates = vertexs[edges[edgeNotPointingToFixedVertex].vertexEnd].Coordinates;
			var a = seeSide(startVertexCoordinates, endVertexCoordinates, fixedCoordinates);
			var b = seeSide(startVertexCoordinates, endVertexCoordinates, newCoordinates);
			var c = seeSide(fixedCoordinates, newCoordinates, startVertexCoordinates);
			var d = seeSide(fixedCoordinates, newCoordinates, endVertexCoordinates);
			if (a * b < 0 && c * d < 0){
				var nextFaceToSearch = getOtherFaceFromEdge(edgeNotPointingToFixedVertex, faceAroundFixedVertex);
				var foundFace =  recursivelySearchByStabbingLine(edgeNotPointingToFixedVertex, nextFaceToSearch, newCoordinates);
				addPointToFace(nPI, foundFace);
				countT++;
			}
		}
	}
}

function recursivelySearchByStabbingLine(edgeEnteredBy, fI, newCoordinates){
	var verticesIndicesAroundFace = getVIF(fI);
	var verticesCoordinatesAroundFace = []
	for(var j = 0; j < verticesIndicesAroundFace.length; ++j){
		verticesCoordinatesAroundFace[j] = vertexs[verticesIndicesAroundFace[j]].Coordinates;
	}
	if (isInsideTriangle(verticesCoordinatesAroundFace, newCoordinates)) return fI
	var edgesAroundFace = getEEF(fI);
	var edgesToCheckWhetherTheyCrossWithStabbingLine = edgesAroundFace.filter(edge => edge !== edgeEnteredBy);
	for (var i = 0; i < edgesToCheckWhetherTheyCrossWithStabbingLine.length; i++){
		var edgeToCheck = edgesToCheckWhetherTheyCrossWithStabbingLine[i];
		var startVertexCoordinates = vertexs[edges[edgeToCheck].vertexStart].Coordinates;;
		var endVertexCoordinates = vertexs[edges[edgeToCheck].vertexEnd].Coordinates;
		var a = seeSide(startVertexCoordinates, endVertexCoordinates, fixedCoordinates);
		var b = seeSide(startVertexCoordinates, endVertexCoordinates, newCoordinates);
		var c = seeSide(fixedCoordinates, newCoordinates, startVertexCoordinates);
		var d = seeSide(fixedCoordinates, newCoordinates, endVertexCoordinates);

		if (a * b < 0 && c * d < 0){
			var nextFaceToSearch = getOtherFaceFromEdge(edgeToCheck, fI);
			return recursivelySearchByStabbingLine(edgeToCheck, nextFaceToSearch, newCoordinates);
		}
	}
}


function addPointToFace(nPI, fI){//
	var countEBeforeAdding = countE;
	var countFBeforeAdding = countF;
	var edgesEnclosingFace = getEEF(fI);
	var aux = getSVFERL(edgesEnclosingFace[0], fI);
	createEdge(countEBeforeAdding,nPI,aux,countFBeforeAdding,fI,countEBeforeAdding + 1,edgesEnclosingFace[2]);
	newOne(edgesEnclosingFace[0], fI, countFBeforeAdding);
	replaceEdge(edgesEnclosingFace[0], edgesEnclosingFace[2], countEBeforeAdding);
	aux = getSVFERL(edgesEnclosingFace[1], fI);
	createEdge(countEBeforeAdding + 1,nPI,aux,countFBeforeAdding + 1,countFBeforeAdding,countEBeforeAdding + 2,edgesEnclosingFace[0]);
	newOne(edgesEnclosingFace[1], fI, countFBeforeAdding + 1);
	replaceEdge(edgesEnclosingFace[1], edgesEnclosingFace[0], countEBeforeAdding + 1);
	aux = getSVFERL(edgesEnclosingFace[2], fI);
	createEdge(countEBeforeAdding + 2,nPI,aux,fI,countFBeforeAdding + 1,countEBeforeAdding,edgesEnclosingFace[1]);
	replaceEdge(edgesEnclosingFace[2], edgesEnclosingFace[1], countEBeforeAdding + 2);
	vertexs[nPI].eI = countEBeforeAdding; 
	faces[fI] = countEBeforeAdding;
	faces[countFBeforeAdding] = countEBeforeAdding + 1; 
	faces[countFBeforeAdding] = countEBeforeAdding + 1; 
	faces[countFBeforeAdding + 1] = countEBeforeAdding + 2;
	countE += 3;
	countF += 2;
}

function computeEnclosingTriangle(pointsE){
    var pointmaxY = pointsE[0];
    var pointminY = pointsE[0];
    for(var i = 0; i < pointsE.length; i++){
        if (pointsE[i].y > pointmaxY.y){
            pointmaxY = pointsE[i]
        }
        if (pointsE[i].y < pointminY.y){
            pointminY = pointsE[i]
        }
    }
    var toAddInY = pointmaxY.y - pointminY.y;
    if (toAddInY === 0) {
        toAddInY = 1;
    }
    var firstVertexEnclosingTriangle = {
        'x': pointmaxY.x,
        'y': pointmaxY.y + toAddInY,
        'z': pointmaxY.z
    };
    var extremePointsVisibleFromFirstVertex = computeExtremePoints(firstVertexEnclosingTriangle, pointsE);
    var mostLeftPoint = extremePointsVisibleFromFirstVertex[0];
    var mostRightPoint = extremePointsVisibleFromFirstVertex[1];
    var secondVertexEnclosingTriangle, thirdVertexEnclosingTriangle;
    if(mostLeftPoint.x === firstVertexEnclosingTriangle.x){
        secondVertexEnclosingTriangle = {
            "x": mostLeftPoint.x + 1, "y": pointminY.y - 1, 'z': 0
        }
    }
    else {
        var leftTangentSlope = (firstVertexEnclosingTriangle.y - mostLeftPoint.y) / (firstVertexEnclosingTriangle.x - mostLeftPoint.x);
        var leftBCoefficient = firstVertexEnclosingTriangle.y - leftTangentSlope * pointmaxY.x;
        var leftTangentXCoordinate = (pointminY.y - leftBCoefficient)/ leftTangentSlope;
        secondVertexEnclosingTriangle = {
            'x': leftTangentXCoordinate + 1, 'y': pointminY.y + 0.3 * leftTangentSlope, 'z': 0
        }
    }
    if(mostRightPoint.x === firstVertexEnclosingTriangle.x){
        thirdVertexEnclosingTriangle = {
            "x": mostRightPoint.x - 1, "y": pointminY.y - 1, "z": 0
        };
    }
    else {
        var rightTangentSlope = (firstVertexEnclosingTriangle.y - mostRightPoint.y) / (firstVertexEnclosingTriangle.x - mostRightPoint.x);
        var rightbCoefficient = firstVertexEnclosingTriangle.y - rightTangentSlope * pointmaxY.x;
        var rightTangentXCoordinate = (pointminY.y - rightbCoefficient)/ rightTangentSlope;
        thirdVertexEnclosingTriangle = {'x': rightTangentXCoordinate - 1, 'y': pointminY.y - 0.3 * rightTangentSlope, 'z': 0}
    }
    var t = [
        firstVertexEnclosingTriangle,
        thirdVertexEnclosingTriangle,
        secondVertexEnclosingTriangle,
    ];
    return t;
}

function seeSide(s1, s2, c){//
	var a = s1;
	var b = s2;
	var d1 = b.x - a.x;
	var d2 = c.y - a.y;
	var d3 = c.x - a.x;
	var d4 = b.y - a.y;
	return (d1)*(d2)-(d3)*(d4);
	//return > 0 left turn, < 0 right turn, == 0 collinear
}

function getTriangleInCounterClockWiseOrder(triangle){//
    var ot = seeSide(triangle[0], triangle[1], triangle[2]);
    if (ot < 0){
        return [triangle[0], triangle[2], triangle[1]]
    }
    return triangle;
}

function isInsideTriangleByDeterminants(ot1, ot2, ot3){//
    return ot1 > 0 && ot2 > 0 && ot3 > 0
}

function isInsideTriangle(triangle, Coordinates){
    var triangleInCounterClockwiseOrder = getTriangleInCounterClockWiseOrder(triangle);
    var vertex1 = triangleInCounterClockwiseOrder[0];
    var vertex2 = triangleInCounterClockwiseOrder[1];
    var vertex3 = triangleInCounterClockwiseOrder[2];
    var ot1 = seeSide(vertex1, vertex2, Coordinates);
    var ot2 = seeSide(vertex2, vertex3, Coordinates);
    var ot3 = seeSide(vertex3, vertex1, Coordinates);
    return ot1 > 0 && ot2 > 0 && ot3 > 0;
}

function isOnBoundaryOfTriangleByDeterminants(ot1, ot2, ot3){
    var sumOfSignFunctionOfseeSide = Math.sign(ot1) + Math.sign(ot2) + Math.sign(ot3);
    return sumOfSignFunctionOfseeSide === 2;
}

function isOnBoundaryOfTriangle(triangle, point){
    var triangleInCounterClockwiseOrder = getTriangleInCounterClockWiseOrder(triangle);
    var vertex1 = triangleInCounterClockwiseOrder[0];
    var vertex2 = triangleInCounterClockwiseOrder[1];
    var vertex3 = triangleInCounterClockwiseOrder[2];
    var ot1 = seeSide(vertex1, vertex2, point);
    var ot2 = seeSide(vertex2, vertex3, point);
    var ot3 = seeSide(vertex3, vertex1, point);
    var sumOfSignFunctionOfseeSide = Math.sign(ot1) + Math.sign(ot2) + Math.sign(ot3);
    return sumOfSignFunctionOfseeSide === 2;
}

function isAVertexOfTriangle(ot1, ot2, ot3){
    var sumOfSignFunctionOfseeSide = Math.sign(ot1) + Math.sign(ot2) + Math.sign(ot3);
    return sumOfSignFunctionOfseeSide === 1;
}

function computeExtremePoints(pointOfView, points){
    var left = points[0];
    var right = points[0];

    for (var i = 1; i < points.length; i++){
        var currentPoint = points[i];
        if (seeSide(pointOfView, left, currentPoint) > 0){
            left = currentPoint;
        }
        else if(seeSide(pointOfView, right, currentPoint) < 0){
            right = currentPoint;
        }
    }
    return [left, right]
}

function computeEnclosingTriangle(pointsE){
    var pointmaxY = pointsE[0];
    var pointminY = pointsE[0];
    for(var i = 0; i < pointsE.length; i++){
        if (pointsE[i].y > pointmaxY.y){
            pointmaxY = pointsE[i]
        }
        if (pointsE[i].y < pointminY.y){
            pointminY = pointsE[i]
        }
    }
    var toAddInY = pointmaxY.y - pointminY.y;
    if (toAddInY === 0) {
        toAddInY = 1;
    }
    var firstVertexEnclosingTriangle = {
        'x': pointmaxY.x,
        'y': pointmaxY.y + toAddInY,
        'z': pointmaxY.z
    };
    var extremePointsVisibleFromFirstVertex = computeExtremePoints(firstVertexEnclosingTriangle, pointsE);
    var mostLeftPoint = extremePointsVisibleFromFirstVertex[0];
    var mostRightPoint = extremePointsVisibleFromFirstVertex[1];
    var secondVertexEnclosingTriangle, thirdVertexEnclosingTriangle;
    if(mostLeftPoint.x === firstVertexEnclosingTriangle.x){
        secondVertexEnclosingTriangle = {
            "x": mostLeftPoint.x + 1, "y": pointminY.y - 1, 'z': 0
        }
    }
    else {
        var leftTangentSlope = (firstVertexEnclosingTriangle.y - mostLeftPoint.y) / (firstVertexEnclosingTriangle.x - mostLeftPoint.x);
        var leftBCoefficient = firstVertexEnclosingTriangle.y - leftTangentSlope * pointmaxY.x;
        var leftTangentXCoordinate = (pointminY.y - leftBCoefficient)/ leftTangentSlope;
        secondVertexEnclosingTriangle = {
            'x': leftTangentXCoordinate + 1, 'y': pointminY.y + 0.3 * leftTangentSlope, 'z': 0
        }
    }
    if(mostRightPoint.x === firstVertexEnclosingTriangle.x){
        thirdVertexEnclosingTriangle = {
            "x": mostRightPoint.x - 1, "y": pointminY.y - 1, "z": 0
        };
    }
    else {
        var rightTangentSlope = (firstVertexEnclosingTriangle.y - mostRightPoint.y) / (firstVertexEnclosingTriangle.x - mostRightPoint.x);
        var rightbCoefficient = firstVertexEnclosingTriangle.y - rightTangentSlope * pointmaxY.x;
        var rightTangentXCoordinate = (pointminY.y - rightbCoefficient)/ rightTangentSlope;
        thirdVertexEnclosingTriangle = {'x': rightTangentXCoordinate - 1, 'y': pointminY.y - 0.3 * rightTangentSlope, 'z': 0}
    }

    var t = [
        firstVertexEnclosingTriangle,
        thirdVertexEnclosingTriangle,
        secondVertexEnclosingTriangle,
    ];
    return t;
}
