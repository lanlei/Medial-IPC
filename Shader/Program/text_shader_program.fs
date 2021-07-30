#ifdef GL_ES
// set default precision to medium
precision mediump int;
precision mediump float;
#endif

uniform sampler2D text;
uniform vec3 textColor;
varying vec2 texcoord;

void main()
{	
	vec4 sampled = vec4(1.0, 1.0, 1.0, texture2D(text, texcoord).r);
	gl_FragColor = vec4(textColor, 1.0)* sampled;
}