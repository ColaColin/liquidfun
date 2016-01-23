/*
* Copyright (c) 2006-2009 Erin Catto http://www.box2d.org
*
* This software is provided 'as-is', without any express or implied
* warranty.  In no event will the authors be held liable for any damages
* arising from the use of this software.
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely, subject to the following restrictions:
* 1. The origin of this software must not be misrepresented; you must not
* claim that you wrote the original software. If you use this software
* in a product, an acknowledgment in the product documentation would be
* appreciated but is not required.
* 2. Altered source versions must be plainly marked as such, and must not be
* misrepresented as being the original software.
* 3. This notice may not be removed or altered from any source distribution.
*/

#include <Box2D/Collision/b2Collision.h>
#include <Box2D/Collision/Shapes/b2PolygonShape.h>
#include <stdio.h>

inline static bool isExtensionConfiguredAsEnd(int32 phantomFlags, int32 index)
{
	return phantomFlags & (1 << (b2_maxPolygonVertices * 2 + index * 2 + 1));
}

// otherCenter is expected to be in the same coordinate system as the vertices of the polygon on which the test is run
inline static bool isExtendedVertex(int32 phantomFlags, int32 index, const b2Vec2* vertices, int32 vCount, const b2Vec2& otherCenter)
{
	bool extended = phantomFlags & (1 << (b2_maxPolygonVertices * 2 + index * 2)); // TODO: if all vertices of a polygon are extended, consider none of them extended
	if (extended)
	{
		bool isEnd = isExtensionConfiguredAsEnd(phantomFlags, index);
		int32 startIndex;
		int32 endIndex;
		if (isEnd)
		{
			startIndex = index == 0 ? vCount-1 : index-1;
			endIndex = index;
		}
		else
		{
			startIndex = index;
			endIndex = (index + 1) % vCount;
		}
		b2Vec2 startVec = vertices[startIndex];
		b2Vec2 endVec = vertices[endIndex];
		return b2PointLineSide(startVec, endVec, otherCenter) < 0;
	}
	else
	{
		return false;
	}
}

inline static b2Vec2 readExtendedVertex(const b2Vec2* vertices, const b2Vec2* normals, int32 count, int32 phantomFlags, int32 index)
{
	b2Vec2 extension;
	if (isExtensionConfiguredAsEnd(phantomFlags, index))
	{
		b2Vec2 normal = normals[(index == 0 ? count-1 : index-1)];
		extension = b2Vec2(-normal.y, normal.x);
	}
	else
	{
		b2Vec2 normal = normals[index];
		extension = b2Vec2(normal.y, -normal.x);
	}
	return vertices[index] + (extension * b2_phantomSlop);
}

// Find the max separation between poly1 and poly2 using edge normals from poly1.
// returns b2_maxFloat if the separation is larger than minSepForSkip
static float32 b2FindMaxSeparation(int32* edgeIndex,
								 const b2PolygonShape* poly1, const b2Transform& xf1,
								 const b2PolygonShape* poly2, const b2Transform& xf2, float32 minSepForSkip, const b2Vec2& poly1Center)
{
	int32 count1 = poly1->m_count;
	int32 count2 = poly2->m_count;
	const b2Vec2* n1s = poly1->m_normals;
	const b2Vec2* n2s = poly2->m_normals;
	int32 p2PhantomFlags = poly2->m_phantomEdges;
	const b2Vec2* v1s = poly1->m_vertices;
	const b2Vec2* v2s = poly2->m_vertices;
	b2Transform xf = b2MulT(xf2, xf1);

	int32 bestIndex = 0;
	float32 trueSeparation = -b2_maxFloat;
	float32 phantomSeparation = -b2_maxFloat;
	
	for (int32 i = 0; i < count1; ++i)
	{
		// Get poly1 normal in frame2.
		b2Vec2 n = b2Mul(xf.q, n1s[i]);
		b2Vec2 v1 = b2Mul(xf, v1s[i]);

		// Find deepest point for normal i.
		float32 si = b2_maxFloat;
		float32 pi = b2_maxFloat;
		for (int32 j = 0; j < count2; ++j)
		{
			float32 sij = b2Dot(n, v2s[j] - v1);
			if (sij < si)
			{
				si = sij;
			}
			
			float32 pij;
			if (isExtendedVertex(p2PhantomFlags, j, v2s, count2, poly1Center))
			{	
				pij = b2Dot(n, readExtendedVertex(v2s, n2s, count2, p2PhantomFlags, j) - v1);
			}
			else
			{
				pij = sij;
			}
			
			if (pij < pi)
			{
				pi = pij;
			}
		}
		
		// no collision will happen anyway, skip the work
		if (si >= minSepForSkip) {
			return b2_maxFloat;
		}
		
		if (si > trueSeparation)
		{
			trueSeparation = si;
		}
		
		if (pi > phantomSeparation)
		{
			phantomSeparation = pi;
			bestIndex = i;
		}
	}
	
	*edgeIndex = bestIndex;
	return phantomSeparation;
}

static void b2FindIncidentEdge(b2ClipVertex c[2],
							 const b2PolygonShape* poly1, const b2Transform& xf1, int32 edge1,
							 const b2PolygonShape* poly2, const b2Transform& xf2, const b2Vec2& center1In2)
{
	const b2Vec2* normals1 = poly1->m_normals;
	int32 p2PhantomFlags = poly2->m_phantomEdges;
	
	int32 count2 = poly2->m_count;
	const b2Vec2* vertices2 = poly2->m_vertices;
	const b2Vec2* normals2 = poly2->m_normals;

	b2Assert(0 <= edge1 && edge1 < poly1->m_count);

	// Get the normal of the reference edge in poly2's frame.
	b2Vec2 normal1 = b2MulT(xf2.q, b2Mul(xf1.q, normals1[edge1]));

	// Find the incident edge on poly2.
	int32 index = 0;
	float32 minDot = b2_maxFloat;
	for (int32 i = 0; i < count2; ++i)
	{
		float32 dot = b2Dot(normal1, normals2[i]);
		if (dot < minDot)
		{
			minDot = dot;
			index = i;
		}
	}

	// Build the clip vertices for the incident edge.
	int32 i1 = index;
	int32 i2 = i1 + 1 < count2 ? i1 + 1 : 0;

	b2Vec2 vi1;
	b2Vec2 vi2;
	if (isExtendedVertex(p2PhantomFlags, i1, vertices2, count2, center1In2))
	{
		vi1 = readExtendedVertex(vertices2, normals2, count2, p2PhantomFlags, i1);
	}
	else
	{
		vi1 = vertices2[i1];
	}
	
	if (isExtendedVertex(p2PhantomFlags, i2, vertices2, count2, center1In2))
	{
		vi2 = readExtendedVertex(vertices2, normals2, count2, p2PhantomFlags, i2);
	}
	else
	{
		vi2 = vertices2[i2];
	}
	
	c[0].v = b2Mul(xf2, vi1);
	c[0].id.cf.indexA = (uint8)edge1;
	c[0].id.cf.indexB = (uint8)i1;
	c[0].id.cf.typeA = b2ContactFeature::e_face;
	c[0].id.cf.typeB = b2ContactFeature::e_vertex;

	c[1].v = b2Mul(xf2, vi2);
	c[1].id.cf.indexA = (uint8)edge1;
	c[1].id.cf.indexB = (uint8)i2;
	c[1].id.cf.typeA = b2ContactFeature::e_face;
	c[1].id.cf.typeB = b2ContactFeature::e_vertex;
}

inline static bool shouldSkipAsPhantom(float32 sep, const b2PolygonShape* p1, int32 edge, const b2Vec2& otherCenter)
{
	if (sep >= -b2_phantomSlop)
	{
		int32 phantomEdges = p1->m_phantomEdges;
		bool vertex0IsPhantom = phantomEdges & (1 << (edge * 2));
		bool vertex1IsPhantom = phantomEdges & (1 << (edge * 2 + 1));
		
		if (vertex0IsPhantom && vertex1IsPhantom) 
		{
			return true;
		}
		if (vertex0IsPhantom || vertex1IsPhantom) 
		{
			const b2Vec2* vertices1 = p1->m_vertices;
			if (vertex0IsPhantom)
			{
				int32 preVI = (edge == 0 ? p1->m_count : edge) - 1;
				b2Vec2 preV = vertices1[preVI];
				b2Vec2 v0 = vertices1[edge];
				if (b2PointLineSide(v0, preV, otherCenter) > 0)
				{
					return true;
				}
			}
			
			if (vertex1IsPhantom)
			{
				b2Vec2 v0 = vertices1[(edge + 1) % p1->m_count];
				b2Vec2 postV = vertices1[(edge + 2) % p1->m_count];
				if (b2PointLineSide(v0, postV, otherCenter) < 0)
				{
					return true;
				}
			}
		}
	}
	return false;
}

// Find edge normal of max separation on A - return if separating axis is found
// Find edge normal of max separation on B - return if separation axis is found
// Choose reference edge as min(minA, minB)
// Find incident edge
// Clip

// The normal points from 1 to 2
void b2CollidePolygons(b2Manifold* manifold,
					  const b2PolygonShape* polyA, const b2Transform& xfA,
					  const b2PolygonShape* polyB, const b2Transform& xfB)
{
	manifold->pointCount = 0;
	float32 totalRadius = polyA->m_radius + polyB->m_radius;

	b2Vec2 centerAInB = polyA->m_centroid;
	if (polyB->m_phantomEdges != 0)
	{
		centerAInB = b2MulT(xfB, b2Mul(xfA, centerAInB));
	} // else centerAInB will never be used, skip the work
	
	int32 edgeA = 0;
	float32 separationA = b2FindMaxSeparation(&edgeA, polyA, xfA, polyB, xfB, totalRadius, centerAInB);
	if (separationA > totalRadius)
		return;

		
	b2Vec2 centerBInA = polyB->m_centroid;
	if (polyA->m_phantomEdges != 0)
	{
		centerBInA = b2MulT(xfA, b2Mul(xfB, centerBInA));
	} // else centerBInA will never be used, skip the work
	
	int32 edgeB = 0;
	float32 separationB = b2FindMaxSeparation(&edgeB, polyB, xfB, polyA, xfA, totalRadius, centerBInA);
	if (separationB > totalRadius)
		return;

	const b2PolygonShape* poly1;	// reference polygon
	const b2PolygonShape* poly2;	// incident polygon
	b2Transform xf1, xf2;
	int32 edge1;					// reference edge
	uint8 flip;
	const float32 k_tol = 0.1f * b2_linearSlop;
	b2Vec2 center1In2;
	b2Vec2 center2In1;
	
	float32 sep1;
	float32 sep2;
	if (separationB > separationA + k_tol)
	{
		poly1 = polyB;
		poly2 = polyA;
		center1In2 = centerBInA;
		center2In1 = centerAInB;
		sep1 = separationB;
		sep2 = separationA;
		xf1 = xfB;
		xf2 = xfA;
		edge1 = edgeB;
		manifold->type = b2Manifold::e_faceB;
		flip = 1;
	}
	else
	{
		poly1 = polyA;
		poly2 = polyB;
		center1In2 = centerAInB;
		center2In1 = centerBInA;
		sep1 = separationA;
		sep2 = separationB;
		xf1 = xfA;
		xf2 = xfB;
		edge1 = edgeA;
		manifold->type = b2Manifold::e_faceA;
		flip = 0;
	}
	
	// this allows to construct bodies from many tightly packed fixtures and ignore collisions on inner edges
	// this prevents many (not all sadly) phantom collisions when objects made from many tiles slide along each other
	// however if objects are completely inside of each other I'd like the behavior not to change, so only ignore
	// phantom collisions that are not too deep, this way polygons can slide on the skin, but will collide when
	// pressed harder into each other. Additionally in case of perfectly straight connections vertices may be extended for 
	// the collision response, see readExtendedVertex 
	if (shouldSkipAsPhantom(sep1, poly1, edge1, center2In1)) 
	{
		return;
	}
	
	b2ClipVertex incidentEdge[2];
	b2FindIncidentEdge(incidentEdge, poly1, xf1, edge1, poly2, xf2, center1In2);
	
	if (shouldSkipAsPhantom(sep2, poly2, incidentEdge[0].id.cf.indexB, center1In2))
	{
		return;
	}

	int32 count1 = poly1->m_count;
	const b2Vec2* vertices1 = poly1->m_vertices;
	const b2Vec2* normals1 = poly1->m_normals;
	
	int32 iv1 = edge1;
	int32 iv2 = edge1 + 1 < count1 ? edge1 + 1 : 0;
	
	int32 phantomFlags = poly1->m_phantomEdges;
	
	b2Vec2 v11;
	b2Vec2 v12;
	if (isExtendedVertex(phantomFlags, iv1, vertices1, count1, center2In1))
	{
		v11 = readExtendedVertex(vertices1, normals1, count1, phantomFlags, iv1);
	}
	else
	{
		v11 = vertices1[iv1];
	}
	
	if (isExtendedVertex(phantomFlags, iv2, vertices1, count1, center2In1))
	{
		v12 = readExtendedVertex(vertices1, normals1, count1, phantomFlags, iv2);
	}
	else
	{
		v12 = vertices1[iv2];
	}

	b2Vec2 localNormal = normals1[edge1];
	b2Vec2 localTangent = b2Vec2(-localNormal.y, localNormal.x);
	
	b2Vec2 planePoint = 0.5f * (v11 + v12);

	b2Vec2 tangent = b2Mul(xf1.q, localTangent);
	b2Vec2 normal = b2Cross(tangent, 1.0f);
	
	v11 = b2Mul(xf1, v11);
	v12 = b2Mul(xf1, v12);

	// Face offset.
	float32 frontOffset = b2Dot(normal, v11);

	// Side offsets, extended by polytope skin thickness.
	float32 sideOffset1 = -b2Dot(tangent, v11) + totalRadius;
	float32 sideOffset2 = b2Dot(tangent, v12) + totalRadius;

	// Clip incident edge against extruded edge1 side edges.
	b2ClipVertex clipPoints1[2];
	b2ClipVertex clipPoints2[2];
	int np;

	// Clip to box side 1
	np = b2ClipSegmentToLine(clipPoints1, incidentEdge, -tangent, sideOffset1, iv1);

	if (np < 2)
		return;

	// Clip to negative box side 1
	np = b2ClipSegmentToLine(clipPoints2, clipPoints1,  tangent, sideOffset2, iv2);

	if (np < 2)
	{
		return;
	}
	
	// Now clipPoints2 contains the clipped points.
	manifold->localNormal = localNormal;
	manifold->localPoint = planePoint;

	int32 pointCount = 0;
	for (int32 i = 0; i < b2_maxManifoldPoints; ++i)
	{
		float32 separation = b2Dot(normal, clipPoints2[i].v) - frontOffset;

		if (separation <= totalRadius)
		{
			b2ManifoldPoint* cp = manifold->points + pointCount;
			cp->localPoint = b2MulT(xf2, clipPoints2[i].v);
			cp->id = clipPoints2[i].id;
			if (flip)
			{
				// Swap features
				b2ContactFeature cf = cp->id.cf;
				cp->id.cf.indexA = cf.indexB;
				cp->id.cf.indexB = cf.indexA;
				cp->id.cf.typeA = cf.typeB;
				cp->id.cf.typeB = cf.typeA;
			}
			++pointCount;
		}
	}

	manifold->pointCount = pointCount;
}
