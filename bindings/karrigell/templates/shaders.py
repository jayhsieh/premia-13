from HTMLTags import SCRIPT

def shader_fs():
   body = """
      #ifdef GL_ES
      precision highp float;
      #endif
      
      varying vec4 vColor;
      varying vec3 vLightWeighting;
      
      void main(void)
      {
         gl_FragColor = vec4(vColor.rgb * vLightWeighting, vColor.a);
      }   
   """
   return SCRIPT(body, id="shader-fs", type="x-shader/x-fragment")

def shader_vs():
   body = """
      attribute vec3 aVertexPosition;
      attribute vec3 aVertexNormal;
      attribute vec4 aVertexColor;
      
      uniform mat4 uMVMatrix;
      uniform mat4 uPMatrix;
      uniform mat3 uNMatrix;
      varying vec4 vColor;
      
      uniform vec3 uAmbientColor;
      uniform vec3 uLightingDirection;
      uniform vec3 uDirectionalColor;
      varying vec3 vLightWeighting;
      
      void main(void)
      {
          gl_Position = uPMatrix * uMVMatrix * vec4(aVertexPosition, 1.0);
          
          vec3 transformedNormal = uNMatrix * aVertexNormal;
          float directionalLightWeighting = max(dot(transformedNormal, uLightingDirection), 0.0);
          vLightWeighting = uAmbientColor + uDirectionalColor * directionalLightWeighting; 

          vColor = aVertexColor;
      }
   """
   return SCRIPT(body, id="shader-vs", type="x-shader/x-vertex")

def axes_shader_fs():
   body = """
      precision mediump float;
      varying vec4 vColor;
      
      void main(void)
      {
         gl_FragColor = vColor;
      }
   """
   return SCRIPT(body, id="axes-shader-fs", type="x-shader/x-fragment")

def axes_shader_vs():
   body = """
      attribute vec3 aVertexPosition;
      attribute vec4 aVertexColor;
      uniform mat4 uMVMatrix;
      uniform mat4 uPMatrix;
      varying vec4 vColor;
      uniform vec3 uAxesColour;
      
      void main(void)
      {
         gl_Position = uPMatrix * uMVMatrix * vec4(aVertexPosition, 1.0);
         vColor =  vec4(uAxesColour, 1.0);
      } 
   """
   return SCRIPT(body, id="axes-shader-vs", type="x-shader/x-vertex")        

def  texture_shader_fs():
   body = """
      #ifdef GL_ES
      precision highp float;
      #endif
      
      varying vec2 vTextureCoord;
      
      uniform sampler2D uSampler;
      
      void main(void)
      {
          gl_FragColor = texture2D(uSampler, vTextureCoord);
      }
   """
   return SCRIPT(body, id="texture-shader-fs", type="x-shader/x-fragment")

def texture_shader_vs():
   body = """
      attribute vec3 aVertexPosition;
      
      attribute vec2 aTextureCoord;
      varying vec2 vTextureCoord;
      
      uniform mat4 uMVMatrix;
      uniform mat4 uPMatrix;
      
      void main(void)
      {
          gl_Position = uPMatrix * uMVMatrix * vec4(aVertexPosition, 1.0);
          vTextureCoord = aTextureCoord; 
      }
   """
   return SCRIPT(body, id="texture-shader-vs", type="x-shader/x-vertex")

print Sum([
      js_include('SurfacePlot.js'),
      js_include('ColourGradient.js'),
      js_include('glMatrix-0.9.5.min.js'),
      js_include('webgl-utils.js'),
      shader_fs(),
      shader_vs(),
      axes_shader_fs(),
      axes_shader_vs(),
      texture_shader_fs(),
      texture_shader_vs(),
   ])
