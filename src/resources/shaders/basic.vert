#version 330 core

layout(location = 0) in vec3 aPos;
layout(location = 1) in float aEdgeT; // для пунктира
uniform mat4 uMVP;
uniform float uPointSize; // размер вершин

out float vEdgeT;

void main() {
  gl_Position = uMVP * vec4(aPos, 1.0);
  vEdgeT = aEdgeT; //VBO
  gl_PointSize = uPointSize;
}
