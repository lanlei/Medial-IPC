#version 330 core
out vec4 FragColor;

in vec3 fragPos;
in vec2 texcoord;
in vec4 lightSpacePos;

uniform vec3 viewPos;

struct DirectionalLight
{
   vec3 position;
   vec3 direction;
   vec4 ambientColor;
   vec4 diffuseColor;
   vec4 specularColor;
};

struct PointLight
{
   vec3 position;
   vec4 ambientColor;
   vec4 diffuseColor;
   vec4 specularColor;
   	float constant;
	float linear;
	float quadratic;
};

#define MAX_LIGHTS_NUM 3 

// The point Light0 is main lighting for calculating shadow.
uniform int pointLightNum;
uniform PointLight pointLight[MAX_LIGHTS_NUM]; //

uniform int directionalLightNum;
uniform DirectionalLight directionalLight[MAX_LIGHTS_NUM]; // 

struct Materials
{
	vec4 ambient;
	vec4 diffuse;
	vec4 specular;
	float shinness;

	int useAmbientMap;
	int useDiffuseMap;
	int useSpecularMap;
	int useBumpMap;
};
uniform Materials material;

uniform int enableLineMode;
uniform int enableShadow;
uniform int useShadowMap;

uniform sampler2D shadowMap; // 0
uniform sampler2D ambientMap; // 1
uniform sampler2D diffuseMap; // 2
uniform sampler2D specularMap; // 3
uniform sampler2D bumpMap; // 4

vec4 _ambientColor;
vec4 _diffuseColor;
vec4 _specularColor;
vec4 _bumpColor;
float _shinness;

vec3 _normal;
vec3 _viewDir;

vec3 CalculateDirLight(DirectionalLight light, vec3 normal, vec3 viewDir);
vec3 CalculatePointLight(PointLight light, vec3 pos, vec3 normal, vec3 viewDir);
float ShadowCalculation(vec4 fragPosLightSpace, vec3 normal, vec3 lightDir);

float rand1(vec2 co);
float rand2(vec2 co);

void main()
{
    if(enableLineMode == 1)
	{
		FragColor = vec4(0.1, 0.1, 0.1, 1);
		return;
	}
	
	_viewDir = normalize(viewPos - fragPos);
	
	if(material.useBumpMap == 1)
	{
		_normal = texture2D(bumpMap, texcoord).rgb;
		_normal = normalize(_normal * 2.0 - 1.0);
	}else _normal = normalize(cross(dFdx(fragPos), dFdy(fragPos)));

	_normal = normalize(cross(dFdx(fragPos), dFdy(fragPos)));
	//
	_ambientColor = material.ambient;
	_diffuseColor = material.diffuse;
	_specularColor = material.specular;
	_bumpColor = vec4(0, 0, 0, 1.0);
	_shinness = material.shinness;

	if(material.useAmbientMap == 1)
		_ambientColor = texture2D(ambientMap, texcoord).rgba;
	if(material.useDiffuseMap == 1)
		_diffuseColor = texture2D(diffuseMap, texcoord).rgba;
	if(material.useSpecularMap == 1)
		_specularColor = texture2D(specularMap, texcoord).rgba;
	//
	
	vec3 result = vec3(0, 0, 0);
	for(int i = 0; i < pointLightNum; i++)
		result += CalculatePointLight(pointLight[i], fragPos, _normal, _viewDir);
	for(int i = 0; i < directionalLightNum; i++)
		result += CalculateDirLight(directionalLight[i], _normal, _viewDir);
		
	float shadow = 1.0;
	if(useShadowMap == 1 && enableShadow == 1)
		shadow = ShadowCalculation(lightSpacePos, _normal, normalize(pointLight[0].position - fragPos));

	result *= shadow;
	FragColor = vec4(result, _ambientColor.a);
}


vec3 CalculateDirLight(DirectionalLight light, vec3 normal, vec3 viewDir)
{
	vec3 lightDir = normalize(-light.direction);
	float diff = max(dot(normal, lightDir), 0.0);

	vec3 halfwayDir = normalize(lightDir + viewDir); 
	float spec = pow(max(dot(viewDir, halfwayDir), 0.0), 64);

	vec3 result = vec3(0, 0, 0);
	result += light.ambientColor.rgb * _ambientColor.rgb;
	result += light.diffuseColor.rgb * diff * _diffuseColor.rgb;
	result += light.specularColor.rgb * spec * _specularColor.rgb;

	return result;
}

vec3 CalculatePointLight(PointLight light, vec3 pos, vec3 normal, vec3 viewDir)
{
	vec3 lightDir = normalize(light.position - pos);
	float diff = max(dot(normal, lightDir), 0.0);

	vec3 halfwayDir = normalize(lightDir + viewDir); 
	float spec = pow(max(dot(viewDir, halfwayDir), 0.0), _shinness);

	//Attenuation
	float distance = length(light.position - fragPos);
	float attenuation = 1.0f / (light.constant + light.linear * distance + light.quadratic * (distance * distance));
	vec3 result = vec3(0, 0, 0);
	result += (light.ambientColor.rgb * _ambientColor.rgb) * attenuation;
	result += (light.diffuseColor.rgb * diff * _diffuseColor.rgb) * attenuation;
	result += (light.specularColor.rgb * spec * _specularColor.rgb) * attenuation;
	return result;
}

float ShadowCalculation(vec4 fragPosLightSpace, vec3 normal, vec3 lightDir)
{
	vec4 shadowMapPosition = fragPosLightSpace / fragPosLightSpace.w;
	shadowMapPosition = (shadowMapPosition + 1.0) /2.0;
	float visibility = 0.0;
	for(int i=-2; i<=2; i++)
      for(int j=-2; j<=2; j++)
        {
            vec2 pos=shadowMapPosition.xy+(vec2(i,j) + vec2(rand1(shadowMapPosition.xy), rand2(shadowMapPosition.xy))*1)*0.005;
		   float distanceFromLight =  texture2D(shadowMap, pos).r; 
            if(distanceFromLight  <  shadowMapPosition.z)
                visibility = visibility + 0.8;
             else visibility = visibility +1.0;    
        }
	visibility=visibility/25.0;
    return visibility;
}

float rand1(vec2 co)
{
	return fract(sin(dot(co.xy,vec2(12.9898,78.233))) * 43758.5453);
}

float rand2(vec2 co)
{
	return fract(sin(dot(co.xy,vec2(12.9898,68.233))) * 33058.5453);
}










