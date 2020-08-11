#pragma once

#include <string>
#include <vector>

// [Input]
// V_in
//     (#V x 3) XYZ coordinate for the mesh vertices
// F_in
//     (#F x {3 or 4}) vertex indices into V ([0, #V-1])
// UV_in
//     (#F x {3 or 4}) UV coordinates of vertices of each face
// VC_in
//     (#V x 4) normalized polypaints (a.k.a vertex color) stored as RGBA ([0.0, 1.0])
//     It seems that alpha values are always 0.0
// M_in
//     (#V x 1) normalized mask ([0.0 (masked), 1.0 (not masked)])
// G_in
//     (#F x 1) face groups stored as a set of ramdom unique ids
// [Output]
// V_out
//     (#V x 3) XYZ coordinate for the mesh vertices
// F_out
//     (#F x {3 or 4}) vertex indices into V ([0, #V-1])
// UV_out
//     (#F x {3 or 4}) UV coordinates of vertices of each face
// VC_out
//     (#V x 4) normalized polypaints (a.k.a vertex color) stored as RGBA ([0.0, 1.0])
//     It seems that alpha values are always 0.0
// M_out
//     (#V x 1) normalized mask ([0.0 (masked), 1.0 (not masked)])
// G_out
//     (#F x 1) face groups stored as a set of ramdom unique ids

template <typename Scalar, typename Index>
void repairAs3DBuilder(
	const std::vector<std::vector<Scalar>> &V_in,
	const std::vector<std::vector<Index>> &F_in,
	const std::vector<std::vector<std::pair<Scalar, Scalar>>> &UV_in,
	const std::vector<std::vector<Scalar>> &VC_in,
	const std::vector<Scalar> &M_in,
	const std::vector<Index> &G_in,
	std::vector<std::vector<Scalar>> &V_out,
	std::vector<std::vector<Index>> &F_out,
	std::vector<std::vector<std::pair<Scalar, Scalar>>> &UV_out,
	std::vector<std::vector<Scalar>> &VC_out,
	std::vector<Scalar> &M_out,
	std::vector<Index> &G_out);

template <typename Scalar, typename Index>
void repairAs3DBuilder(
	const std::vector<std::vector<Scalar>> &V_in,
	const std::vector<std::vector<Index>> &F_in,
	std::vector<std::vector<Scalar>> &V_out,
	std::vector<std::vector<Index>> &F_out);

#include "repairAs3DBuilder.cpp"
