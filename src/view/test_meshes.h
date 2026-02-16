#pragma once

#include "../render/opengl/gl_render.h"
namespace s21::test_meshes {

	// Каркасный куб: 8 вершин, 12 рёбер (24 индекса)
	inline GlRender::MeshData MakeWireCube(float half = 0.5f) {
		GlRender::MeshData m;

		const float v[] = {
			-half, -half, -half,  // 0
			 half, -half, -half,  // 1
			 half,  half, -half,  // 2
			-half,  half, -half,  // 3
			-half, -half,  half,  // 4
			 half, -half,  half,  // 5
			 half,  half,  half,  // 6
			-half,  half,  half   // 7
		};
		m.vertices_xyz.assign(v, v + 8 * 3);

		const std::uint32_t e[] = {
			0,1, 1,2, 2,3, 3,0,
			4,5, 5,6, 6,7, 7,4,
			0,4, 1,5, 2,6, 3,7
		};
		m.edge_indices.assign(e, e + 24);

		// Для заливки (если включите fill): 12 треугольников = 36 индексов
		const std::uint32_t t[] = {
			0,1,2, 0,2,3,
			4,6,5, 4,7,6,
			0,3,7, 0,7,4,
			1,5,6, 1,6,2,
			0,4,5, 0,5,1,
			3,2,6, 3,6,7
		};
		m.tri_indices.assign(t, t + 36);

		return m;
	}
	GlRender::MeshData MakeInnerFrameCube(float outer_half = 0.75f, float inner_half = 0.25f);

	GlRender::MeshData MakeTorus(float majorR, float minorR, int majorSeg, int minorSeg);


}

