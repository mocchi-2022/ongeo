// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#include "ONGEO.h"

#include <cstdio>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <algorithm>

bool ONGEO_ReadBinarySTL(ON_BinaryArchive &ba, ON_Mesh &mesh, ON_String &header){
	mesh.Destroy();
	ba.SeekFromStart(0);
	header.ReserveArray(80);
	if (!ba.ReadByte(80, header.Array())) return false;
	int num_triangles;
	if (!ba.ReadInt(&num_triangles)) return false;
	for (int i = 0, i3 = 0; i < num_triangles; ++i, i3 += 3){
		float v[12];
		if (!ba.ReadFloat(12, v)) break;
		mesh.m_FN.Append(ON_3fVector(v));
		mesh.SetVertex(i3, ON_3fPoint(v+3));
		mesh.SetVertex(i3+1, ON_3fPoint(v+6));
		mesh.SetVertex(i3+2, ON_3fPoint(v+9));
		mesh.SetTriangle(i, i3, i3+1, i3+2);
		ba.SeekFromCurrentPosition(2);
	}
	mesh.CombineIdenticalVertices(false, true);
	return true;
}

void ONGEO_ReadLine(ON_BinaryArchive &ba, ON_String &str){
	str.Empty();
	char buf[100];
	for (;;){
		size_t bs = 100;
		size_t cp = ba.CurrentPosition();
		if (!ba.ReadChar(100, buf)){
			bs = ba.CurrentPosition() - cp;
			ba.SeekFromCurrentPosition(-static_cast<int>(bs));
			ba.ReadChar(bs, buf);
		}
		for (size_t i = 0; i < bs; ++i){
			if (buf[i] != '\r' && buf[i] != '\n') continue;

			str.Append(buf, i);
			size_t nl = i + 1;
			if (nl < bs && buf[nl-1] == '\r' && buf[nl] == '\n') ++nl;
			ba.SeekFromCurrentPosition(nl-bs);
			return;
		}
		str.Append(buf, bs);
	}
}

bool ONGEO_ReadTextSTL(ON_BinaryArchive &ba, ON_Mesh &mesh, ON_String &modelname){
	mesh.Destroy();
	ON_String line;
	ONGEO_ReadLine(ba, line);
	modelname.ReserveArray(line.Length() - 6);
	if (std::sscanf(line, "solid %s", modelname.Array()) == 0) return false;
	int linestate = 0;
	ON_3fVector fn;
	ON_3fPoint v[3];
	bool suspicion_error = false;
	for(;!ba.AtEnd();){
		ONGEO_ReadLine(ba, line);
		if (line.IsEmpty()) continue;
		
		if (suspicion_error){
			for (int i = 0; i < line.Length(); ++i){
				if (line[i] < 0 || line[i] >= 256 || !std::isalnum(line[i])) return false;
			}
		}

		switch(linestate){
		case 0:
			if (suspicion_error = (std::sscanf(line, "%*[ \t]facet normal %f %f %f", &fn.x, &fn.y, &fn.z) < 3))
				continue;
			++linestate;
			break;
		case 1:
			if (suspicion_error = (std::strstr(line, "outer loop") == 0)) continue;
			++linestate;
			break;
		case 2:
		case 3:
		case 4:
			if (suspicion_error = (std::sscanf(line, "%*[ \t]vertex %f %f %f", &v[linestate-2].x, &v[linestate-2].y, &v[linestate-2].z) < 3)) continue;
			++linestate;
			break;
		case 5:
			if (suspicion_error = (std::strstr(line, "endloop") == 0)) continue;
			++linestate;
			break;
		case 6:
			if (suspicion_error = (std::strstr(line, "endfacet") == 0)) continue;
			linestate = 0;
			mesh.m_FN.Append(fn);
			mesh.m_V.Append(v[0]);
			mesh.m_V.Append(v[1]);
			mesh.m_V.Append(v[2]);
			int tn = mesh.m_FN.Count()-1;
			mesh.SetTriangle(tn, tn*3, tn*3+1, tn*3+2);
			break;
		}
	}
	mesh.CombineIdenticalVertices();
	return true;
}
