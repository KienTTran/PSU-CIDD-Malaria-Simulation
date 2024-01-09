//
// Created by kient on 6/16/2023.
//

#ifndef MASS_RENDERENTITY_CUH
#define MASS_RENDERENTITY_CUH

#include "GL/glew.h"
#include "Helpers/Shader.h"
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec2.hpp>

namespace GPU{
    class RenderEntity;
}

class Model;

typedef struct {
    GLuint count;
    GLuint primCount;
    GLuint firstIndex;
    GLuint baseVertex;
    GLuint baseInstance;
} DrawElementsIndirectCommand;

class GPU::RenderEntity {
public:
    RenderEntity(Model *model = nullptr);
    ~RenderEntity();

    void init_entity();
    void init_render_entity(int window_width, int window_height);

public:
    Model* model_;
    GLuint VAO;
    Shader *shader;

    std::vector<glm::vec4> entity_vertices;
    std::vector<glm::vec4> entity_colors;
    glm::mat4 projection;
    glm::mat4 view;
    GLint entity_indices[3] = {0,1,2};
//    GLint entity_indices[6] = {0,1,2,0,1,3};
    GLuint *VBO;
    GLuint EBO;
    GLuint *SSBO;
    GLuint CMD;
    //On DEVICE OGL
    struct cudaGraphicsResource *d_cuda_buffer_model;
    size_t d_ogl_buffer_model_num_bytes; // to get models data from gpu_buffer
    struct cudaGraphicsResource *d_cuda_buffer_color;
    size_t d_ogl_buffer_color_num_bytes;// to get colors data from gpu_buffer
    glm::mat4 *d_ogl_buffer_model_ptr;
    glm::vec4 *d_ogl_buffer_color_ptr;
};

#endif //MASS_RENDERENTITY_CUH
