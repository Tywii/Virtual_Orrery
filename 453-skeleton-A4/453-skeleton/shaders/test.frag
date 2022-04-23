#version 330 core

in vec3 fragPos;
in vec3 fragColor;
in vec3 n;
in vec2 tc;

uniform vec3 light;
uniform sampler2D sampler;

out vec4 color;

void main() {
	vec4 d = texture(sampler, tc);
	vec3 lightDir = normalize(light - fragPos);
    vec3 normal = normalize(n);
    float diff = max(dot(lightDir, normal), 0.0);
    //float diff = max(dot(normal, lightDir), 0.0);

	//vec3 result = diff * d;
	color = 0.5 * d; 
}
