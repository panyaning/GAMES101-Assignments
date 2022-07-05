#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
{
    Vector2D node_distance = (end-start)/(num_nodes+1);

    for (int i= 0; i <num_nodes+2; i++){
        Mass *node = new Mass(start+node_distance*i, node_mass, false);
        masses.push_back(node);
        if(i ==0) continue;
        Spring* spring =new Spring(masses[i-1], masses[i], k);
        springs.push_back(spring);
    }
    for (auto &i : pinned_nodes) {
        masses[i]->pinned = true;
    }
}

void Rope::simulateEuler(float delta_t, Vector2D gravity)
{
    for (auto &s : springs)
    {
        Vector2D v2_ba = s->m2->position - s->m1->position;
        double dis = v2_ba.norm();
        Vector2D dir_ba = v2_ba / dis;
        Vector2D f_ab = s->k
                * dir_ba
                * (dis  - s->rest_length);
        Vector2D f_ba = -f_ab;

        s->m1->forces +=f_ab;
        s->m2->forces +=f_ba;
    }
    for (auto &m : masses)
    {
        if (!m->pinned)
        {
            Vector2D a = m->forces /m->mass + gravity;
            m->velocity += a *delta_t;
            //m->position += m->velocity *delta_t;//显示欧拉
            Vector2D v_new = m->velocity + a * delta_t;
            m->position += v_new * delta_t; //半隐式
            //隐式欧拉不会
        }
        m->forces = Vector2D(0, 0);
    }
}

void Rope::simulateVerlet(float delta_t, Vector2D gravity)
{
    for (auto &s : springs)
    {
        Vector2D v2_ba = s->m2->position - s->m1->position;
        double dis = v2_ba.norm();
        Vector2D dir_ba = v2_ba / dis;

        if (s->m1->pinned && s->m2->pinned) continue; //m1写成m2 找了一小时......
        if ( s->m2->pinned && !s->m1->pinned){//m1 不是固定点 m2是
            s->m1->position += dir_ba *(dis  - s->rest_length) ;
        }else if( s->m1->pinned && !s->m2->pinned){//m1是固定点 m2不是
            s->m2->position += -dir_ba *(dis  - s->rest_length);
        }else {//m1m2不是
            s->m1->position += dir_ba *(dis  - s->rest_length) * 0.5;
            s->m2->position += -dir_ba *(dis  - s->rest_length) * 0.5;
        }
    }
    for (auto &m : masses)
    {
        if (!m->pinned)
        {
            Vector2D a = m->forces /m->mass + gravity;
            Vector2D v_old = m->position;
            double damping = 0.00005;
            //显式 Verlet add damping
            m->position += (1-damping) * (m->position - m->last_position)  + a *delta_t*delta_t;
            //显式 Verlet
            //m->position += m->position - m->last_position + a *delta_t*delta_t;
            m->last_position = v_old;
        }
    }
} 
}
