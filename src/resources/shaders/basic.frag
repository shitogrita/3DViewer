#version 330 core

in float vEdgeT; //параметр для пунктира
uniform vec4 uColor; // текущий цвет

uniform int uUseDash; // тип линии ребра 0 - сплошная, 1 - пунктир
uniform float uDashPeriod; // период пунктира
uniform float uDashFill; // ширина штриха
uniform int uPointMode; // 0 - ничего, 1 - круг, 2 - квадрат
uniform float uPointSoft; // мягкость края точки

out vec4 FragColor;

void main() {
    if (uUseDash == 1) { // пунктир
        float period = max(uDashPeriod, 1e-6);
        float x = fract(vEdgeT / period);
        if (8.0 * x > clamp(uDashFill, 0.0, 1.0)) discard;
    }
    if (uPointMode == 1) {  // круг
        vec2 p = gl_PointCoord - vec2(0.5);
        float r = length(p);
        float edge0 = 0.5 - uPointSoft;
                float edge1 = 0.5;
                float alpha = 1.0 - smoothstep(edge0, edge1, r);
                if (alpha <= 0.0) discard;
                FragColor = vec4(uColor.rgb, uColor.a * alpha);
                return;
    } else if (uPointMode == 2 ) { // квадрат
        vec2  p = abs(gl_PointCoord - vec2(0.5));
                float d = max(p.x, p.y);
                float edge0 = 0.5 - uPointSoft;
                float edge1 = 0.5;
                float alpha = 1.0 - smoothstep(edge0, edge1, d);
                if (alpha <= 0.0) discard;
                FragColor = vec4(uColor.rgb, uColor.a * alpha);
                return;
    }
    FragColor = uColor;
}
