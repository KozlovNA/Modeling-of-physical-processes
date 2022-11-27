#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory>
#include <string>
#include <valarray>

using namespace std;
#include "json.hpp"
#include "json_fwd.hpp"
using json = nlohmann::json;

// parameters of simulation
double duration;
double dt;
int n;
string OUTPATH;
void init(string config){
    std::ifstream i(config, std::ifstream::binary);
        json j;
        i >> j;
    duration = j["StartCondition"]["duration"];
    dt = j["StartCondition"]["dt"];
    n = floor(duration/dt);
    OUTPATH = j["OUTPATH"];
}


//-----------MATERIAL DOT OBJECT-----------//
//--contains: start conditions (x_0, v_0)--//
//------------angular frequency (w)--------//
//------------x,v,energy arrays------------//
class Object {
public:
    double x_0, v_0, w;
    double* x;
    double* v;
    double* energy;
    Object():
        x(new double[n]),
        v(new double[n]),
        energy(new double[n]),
        x_0(10), v_0(0), w(1)
        {   
            std::ifstream i("config.json");
            json j;
            i >> j;
            *x = j["StartCondition"]["x_0"];
            *v = j["StartCondition"]["v_0"];
            w = j["StartCondition"]["w"];
        }

    ~Object()
    {
        delete[] x;
        delete[] v;
        delete[] energy;
    }

    Object& operator=(Object const &src) {

        double *new_x = new double[n];
        for (size_t pos = 0; pos != n; ++pos)
            new_x[pos] = x[pos];
        delete[] x;
        x = new_x;

        double *new_v = new double[n];
        for (size_t pos = 0; pos != n; ++pos)
            new_v[pos] = v[pos];
        delete[] v;
        v = new_v;

        double *new_energy = new double[n];
        for (size_t pos = 0; pos != n; ++pos)
            new_energy[pos] = energy[pos];
        delete[] energy;
        energy = new_energy;

        return *this;
    }

    Object(Object const &src): Object(){
        for (size_t pos = 0; pos != n; pos++)
            x[pos] = src.x[pos];

        for (size_t pos = 0; pos != n; pos++)
            v[pos] = src.v[pos];

        for (size_t pos = 0; pos != n; pos++)
            energy[pos] = src.energy[pos];
    }

    void measureEnergy() {
        for (int i = 0; i < n; i++) {
            energy[i] = v[i] * v[i] / 2 + w * w * x[i] * x[i] / 2;
        }
    }

};


//---Kahan Summation for array of arguments---//
double KahanSum(double* input, int inputsize) {
    double sum = 0.0;
    double c = 0.0;
    for(int i = 1; i < inputsize; i++){
        double y = input[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}


//----------------Kahan Summation for 2 arguments----------------//
//---needs to initiate double c variable for error remembering---//
double KahanSum(double input, double term, double c) {
        double y = term - c;
        double t = input + y;
        c = (t - input) - y;
    return t;
}


//---------------------------------------//
//------NUMERIC CALCULATION METHODS------//
//---------------------------------------//
class HoynsSchemeKahan{
public:
    double* x_1;
    double* v_1;

public:
    HoynsSchemeKahan():
        x_1(new double[n]),
        v_1(new double[n]) {}

    ~HoynsSchemeKahan() {delete[] x_1; delete[] v_1;}

void count(Object &object){
    double c_x = 0.0;
    double c_v = 0.0;
    for (int i = 0; i < n - 1; i++){
        x_1[i + 1] = x_1[i] + dt * v_1[i] - c_x;
        c_x = (x_1[i + 1] - x_1[i]) - (dt * v_1[i] - c_x);
        v_1[i + 1] = v_1[i] + dt * (-object.w * object.w) * x_1[i] - c_v;
        c_v = (v_1[i + 1] - v_1[i]) - (dt * (-object.w * object.w) * x_1[i] - c_v);
    }
    c_x = 0.0;
    c_v = 0.0;
    for (int i = 0; i < n - 1; i++) {
        object.x[i + 1] = object.x[i] + dt/2 * (object.v[i] + v_1[i + 1]) - c_x;
        c_x = (object.x[i + 1] - object.x[i]) - (dt/2 * (object.v[i] + v_1[i + 1]) - c_x);
        object.v[i + 1] = object.v[i] + dt/2 * (-object.w * object.w)*(object.x[i] + x_1[i + 1]) - c_v;
        c_v = (object.v[i + 1] - object.v[i]) - (dt/2 * (-object.w * object.w)*(object.x[i] + x_1[i + 1]) - c_v);
    }
}
     HoynsSchemeKahan& operator=(HoynsSchemeKahan const &src) = delete;
     HoynsSchemeKahan(HoynsSchemeKahan const &src) = delete;
};


class HoynsScheme{
public:
    double* x_1;
    double* v_1;

public:
    HoynsScheme():
            x_1(new double[n]),
            v_1(new double[n]) { }

    ~HoynsScheme() {delete[] x_1; delete[] v_1;}

    void count(Object &object){
        for (int i = 0; i < n - 1; i++){
            x_1[i + 1] = x_1[i] + dt * v_1[i];
            v_1[i + 1] = v_1[i] + dt * (-object.w * object.w) * x_1[i];
        }
        for (int i = 0; i < n - 1; i++) {
            object.x[i + 1] = object.x[i] + dt/2 * (object.v[i] + v_1[i + 1]);
            object.v[i + 1] = object.v[i] + dt/2 * (-object.w * object.w)*(object.x[i] + x_1[i + 1]);
        }
    }

    HoynsScheme& operator=(HoynsScheme const &src) = delete;

    HoynsScheme(HoynsScheme const &src) = delete;
};


class EulersMethod{
public:
    EulersMethod() {}

    void count(Object &object){  
        for (int i = 0; i < n - 1; i++){
            object.x[i + 1] = object.x[i] + dt * object.v [i];
            object.v[i + 1] = object.v[i] + dt * (-object.w*object.w) * object.x[i]; 
        }
    }

    EulersMethod& operator=(EulersMethod const &src) = delete;

    EulersMethod(EulersMethod const &src) = delete;

    ~EulersMethod() {}
};

//---------GENERIC--ALGORYTHMS--------//
class MathPendState{
public:
    array<double, 2> state;
    double w;

    MathPendState() {
        std::ifstream i("config.json", std::ifstream::binary);
        json j;
        i >> j;
        state[0] = j["StartCondition"]["x_0"];
        state[1] = j["StartCondition"]["v_0"];
        w = j["StartCondition"]["w"];
    }

    MathPendState(double x_0, double v_0, double w) {
        state[0] = x_0;
        state[1] = v_0;
        this->w = w;
    }

    MathPendState& operator=(MathPendState const &src){
        state[0] = src.state[0];
        state[1] = src.state[1];
        w = src.w;
        return *this;
    }

    MathPendState(MathPendState const &src){
        state[0] = src.state[0];
        state[1] = src.state[1];
        w = src.w;
    }

    ~MathPendState() {}

    double operator[](int i){
        if (i > 1) throw std::out_of_range("Index is out of range");
        return state[i];
    }

    MathPendState f(MathPendState &state){
        MathPendState result;
        result.state[0] = state.state[1];
        result.state[1] = -w * w * state.state[0];
        return result;
    }

    MathPendState operator+(MathPendState const &src){
        MathPendState result;
        result.state[0] = state[0] + src.state[0];
        result.state[1] = state[1] + src.state[1];
        return result;
    }

    MathPendState operator*(double a){
        MathPendState result;
        result.state[0] = this->state[0] * a;
        result.state[1] = this->state[1] * a;
        return result;
    }

    MathPendState operator/(double a){
        MathPendState result;
        result.state[0] = this->state[0] / a;
        result.state[1] = this->state[1] / a;
        return result;
    }
};


class PhysPendState{
public:
    array<double, 2> state;
    double g;
    double l;
    PhysPendState(){
        std::ifstream i("config.json", std::ifstream::binary);
        json j;
        i >> j;
        state[0] = j["PhysicalPendulum"]["teta_0"];
        state[1] = j["PhysicalPendulum"]["d(teta_0)/dt"];
        g = j["PhysicalPendulum"]["g"];
        l = j["PhysicalPendulum"]["lambda"];
    }
    PhysPendState(double x, double v){
        PhysPendState st;
        g = st.g;
        l = st.l;
        state[0] = x;
        state[1] = v;
    }
    PhysPendState(PhysPendState const &single_state){
        PhysPendState st;
        g = st.g;
        l = st.l;
        state[0] = single_state.state[0];
        state[1] = single_state.state[1];
    }

    PhysPendState& operator=(PhysPendState const &src){
        state[0] = src.state[0];
        state[1] = src.state[1];
        g = src.g;
        l = src.l;
        return *this;
    }

    double operator[](int i){
        if (i > 1) throw std::out_of_range("Index is out of range");
        return state[i];
    }

    PhysPendState f(PhysPendState &state){
        PhysPendState result;
        result.state[0] = state[1];
        result.state[1] = -g/l * sin(state[0]);
        return result;
    }
    
    PhysPendState operator+(PhysPendState const &src){
        PhysPendState result;
        result.state[0] = state[0] + src.state[0];
        result.state[1] = state[1] + src.state[1];
        return result;
    }
    
    PhysPendState operator*(double const &src){
        PhysPendState result;
        result.state[0] = state[0] * src;
        result.state[1] = state[1] * src;
        return result;
    }

    PhysPendState operator/(double const &src){
        PhysPendState result;
        result.state[0] = state[0] / src;
        result.state[1] = state[1] / src;
        return result;
    }

};

class PhysGammaPendState{
public:
    array<double, 2> state;
    double w;
    double y;

    PhysGammaPendState(){
        std::ifstream i("config.json", std::ifstream::binary);
        json j;
        i >> j;
        state[0] = j["PhysicalPendulum2"]["teta_0"];
        state[1] = j["PhysicalPendulum2"]["d(teta_0)/dt"];
        w = j["PhysicalPendulum2"]["w"];
        y = j["PhysicalPendulum2"]["y"];
    }

    PhysGammaPendState(double x, double v){
        PhysGammaPendState st;
        w = st.w;
        y = st.y;
        state[0] = x;
        state[1] = v;
    }
    PhysGammaPendState(PhysGammaPendState const &single_state){
        PhysGammaPendState st;
        w = st.w;
        y = st.y;
        state[0] = single_state.state[0];
        state[1] = single_state.state[1];
    }

    PhysGammaPendState& operator=(PhysGammaPendState const &src){
        state[0] = src.state[0];
        state[1] = src.state[1];
        w = src.w;
        y = src.y;
        return *this;
    }

    double operator[](int i){
        if (i > 1) throw std::out_of_range("Index is out of range");
        return state[i];
    }

    PhysGammaPendState f(PhysGammaPendState &state){
        PhysGammaPendState result;
        result.state[0] = state[1];
        result.state[1] = - w * w * sin(state[0]) - 2 * y * state[1];
        return result;
    }
    
    PhysGammaPendState operator+(PhysGammaPendState const &src){
        PhysGammaPendState result;
        result.state[0] = state[0] + src.state[0];
        result.state[1] = state[1] + src.state[1];
        return result;
    }
    
    PhysGammaPendState operator*(double const &src){
        PhysGammaPendState result;
        result.state[0] = state[0] * src;
        result.state[1] = state[1] * src;
        return result;
    }

    PhysGammaPendState operator/(double const &src){
        PhysGammaPendState result;
        result.state[0] = state[0] / src;
        result.state[1] = state[1] / src;
        return result;
    }
};


class PhysDrivenPendState{
public:
    array<double, 3> state;
    double w;
    double y;
    double F;
    double w_F;
    double sk;

    PhysDrivenPendState(){
        std::ifstream i("config.json", std::ifstream::binary);
        json j;
        i >> j;
        state[0] = j["DrivenForcePendulum"]["teta_0"];
        state[1] = j["DrivenForcePendulum"]["d(teta_0)/dt"];
        state[2] = 0;
        w = j["DrivenForcePendulum"]["w"];
        y = j["DrivenForcePendulum"]["y"];
        F = j["DrivenForcePendulum"]["F"];
        w_F = j["DrivenForcePendulum"]["w_F"];
        sk = j["DrivenForcePendulum"]["skvazh"];
    }

    PhysDrivenPendState(double x, double v, double t){
        PhysDrivenPendState st;
        w = st.w;
        y = st.y;
        F = st.F;
        w_F = st.w_F;
        sk = st.sk;
        state[0] = x;
        state[1] = v;
        state[2] = t;

    }
    PhysDrivenPendState(PhysDrivenPendState const &single_state){
        PhysDrivenPendState st;
        w = st.w;
        y = st.y;
        F = st.F;
        w_F = st.w_F;
        sk = st.sk;
        state[0] = single_state.state[0];
        state[1] = single_state.state[1];
        state[2] = single_state.state[2];
    }

    PhysDrivenPendState& operator=(PhysDrivenPendState const &src){
        state[0] = src.state[0];
        state[1] = src.state[1];
        state[2] = src.state[2];
        w = src.w;
        y = src.y;
        F = src.F;
        w_F = src.w_F;
        sk = src.sk;
        return *this;
    }

    double operator[](int i){
        if (i > 2) throw std::out_of_range("Index is out of range");
        return state[i];
    }

    double fsin(double t){
        return F * sin(w_F * t);
    }

    double fcos(double t){
        return F * cos(w_F * t);
    }

    double meandr(double t){
        if (fmod(t, 2*3.1416/w_F) < 2*3.1416/w_F*(1-sk)) return 1;
        else return 0;
    }

    PhysDrivenPendState f(PhysDrivenPendState &state){
        PhysDrivenPendState result;
        result.state[2] = state[2] + dt;
        result.state[0] = state[1];
        result.state[1] = - w * w * state[0] - 2 * y * state[1] + fcos(state[2]);
        return result;
    }
    
    PhysDrivenPendState operator+(PhysDrivenPendState const &src){
        PhysDrivenPendState result;
        result.state[0] = state[0] + src.state[0];
        result.state[1] = state[1] + src.state[1];
        result.state[2] = src.state[2];
        return result;
    }
    
    PhysDrivenPendState operator*(double const &src){
        PhysDrivenPendState result;
        result.state[0] = state[0] * src;
        result.state[1] = state[1] * src;
        result.state[2] = state[2];
        return result;
    }

    PhysDrivenPendState operator/(double const &src){
        PhysDrivenPendState result;
        result.state[0] = state[0] / src;
        result.state[1] = state[1] / src;
        result.state[2] = state[2];
        return result;
    }
};



template<typename State>
class GenericEuler{
public:
    State step(State &state){
        State new_state;
        new_state = state + state.f(state)*dt;
        return new_state;
    } 
};
template<typename State>
class GenericHoyn{
    public:
    State mezo_step(State &state){
        State mezo_state;
        mezo_state = state + state.f(state)*dt;
        return mezo_state;
    }
    State step(State &state){
        State new_state;
        State mezo_state = mezo_step(state);
        new_state = state + (state.f(state) + state.f(mezo_state))*dt/2;
        return new_state;
    }
};

template<typename State>
class GenericRK{
private:
    State k_1_res;
    State k_2_res;
    State k_3_res;
    State k_4_res;
    
    void k_1(State &state){
        k_1_res = state.f(state)*dt;
    }

    void k_2(State &state){
        State mezo_state;
        mezo_state = state + k_1_res/2.0;
        k_2_res = mezo_state.f(mezo_state)*dt;
    }

    void k_3(State &state){
        State mezo_state;
        mezo_state = state + k_2_res/2.0;
        k_3_res = mezo_state.f(mezo_state)*dt;
    }

    void k_4(State &state){
        State mezo_state;
        mezo_state = state + k_3_res;
        k_4_res = mezo_state.f(mezo_state)*dt;
    }
public:
    State step(State &state){
        k_1(state);
        k_2(state);
        k_3(state);
        k_4(state);
        //cout << "k_1_res" << k_1_res.state[0] << '\n';
        //cout << k_2_res.state[0] << '\n';
        //cout << k_3_res.state[0] << '\n';
        //cout << k_4_res.state[0] << '\n';
        State new_state;
        new_state = state + (k_1_res + k_2_res*2 + k_3_res*2 + k_4_res)/6.0;
        return new_state;
    }

};

class PGP_presise_solution{
public:
    vector<array<double, 2>> solution;
    double w;
    double y;
    double teta_0;
    double d_teta_0;

    PGP_presise_solution(){
        std::ifstream i("config.json", std::ifstream::binary);
        json j;
        i >> j;
        teta_0 = j["PhysicalPendulum2"]["teta_0"];
        d_teta_0 = j["PhysicalPendulum2"]["d(teta_0)/dt"];
        w = j["PhysicalPendulum2"]["w"];
        y = j["PhysicalPendulum2"]["y"];
    }

    void solve(string output){
        for (int i = 0; i < n; i++){
            array<double, 2> temp;
            if (w > y) { 
            temp[0] = teta_0*exp(-y*i*dt)*cos(sqrt(w*w-y*y)*i*dt);
            //temp[1] = teta_0*exp(-y*i*dt)*cos(w*i*dt)
            temp[1] = 0;
            }
            solution.push_back(temp);
        }
        ofstream out;
        out.open(output);
        if (out.is_open())
        {
            for (int i = 0; i < n; i++){
                out <<  dt * i  << ' ' << solution[i][0]<<'\n';
            }
        }
        out.close();
    }
};

class DO_precise_solution{
public:
    vector<array<double, 2>> solution;
    double w;
    double y;
    double F;
    double w_F;
    double sk;
    double phi_0;

    DO_precise_solution(){
        std::ifstream i("config.json", std::ifstream::binary);
        json j;
        i >> j;
        w = j["DrivenForcePendulum"]["w"];
        y = j["DrivenForcePendulum"]["y"];
        F = j["DrivenForcePendulum"]["F"];
        w_F = j["DrivenForcePendulum"]["w_F"];
        sk = j["DrivenForcePendulum"]["skvazh"];
        phi_0 = atan(2*y*w_F/(w_F*w_F - w*w));
    }

    void solve(string output){

        for (int i = 0; i < n; i++){
            array<double, 2> temp;
            if (w > y) {
            temp[0] = F/(sqrt((w*w-w_F*w_F)*(w*w - w_F*w_F)+4*y*y*w_F*w_F))*(cos(w_F*i*dt+phi_0) - exp(-y*i*dt)*cos(phi_0)*cos(sqrt(w*w-y*y)*i*dt) + (w_F*sin(phi_0)-cos(phi_0))/(sqrt(w*w - y*y))*exp(-y*i*dt)*sin(sqrt(w*w-y*y)*i*dt));
            temp[1] = 0; //velocity not done yet
            }
            solution.push_back(temp);
        }
        ofstream out;
        out.open(output);
        if (out.is_open())
        {
            for (int i = 0; i < n; i++){
                out <<  dt * i  << ' ' << solution[i][0]<<'\n';
            }
        }
        out.close();
    }
};

class PGP_Q_finder{
public:
    vector <array<double, 2>> Q;
    double func(){
        
    }
};

template<typename State, typename Method>
void solver(string output)
{
    OUTPATH = output;
    Method method;
    vector<State> states;
    State initstate;
    states.push_back(initstate);
    for (int i = 0; i < n; i++){
        State new_state(method.step(states.back())); 
        states.push_back(new_state);
    }
    ofstream out;
        out.open(OUTPATH);
        if (out.is_open())
        {
            for (int i = 0; i < n; i++){
                out <<  dt * i  << ' ' << states[i].state[0]<< ' ' << states[i].state[1] <<'\n';
            }
        }
        out.close();
}

//--------------------------------------//
class RungeKutta{
public:
    RungeKutta() {}
    ~RungeKutta(){}

    double *f(Object &object, int i){
        static double res[2] = {0, 0};
        res[0] = object.v[i];
        res[1] = -object.x[i]*object.w*object.w;
        return res;
    }

    double *k_1(Object &object, int i){
        static double res[2] = {0, 0};
        res[1] = f(object, i)[1]*dt;
        res[0] = f(object, i)[0]*dt;
        return res;
    }
    double *k_2(Object &object, int i){
        static double res[2] = {0, 0};
        Object object2;
        object2.x[i] = object.x[i] + 1.0/2.0 * k_1(object, i)[0];
        object2.v[i] = object.v[i] + 1.0/2.0 * k_1(object, i)[1];
        res[0] = f(object2, i)[0]*dt;
        res[1] = f(object2, i)[1]*dt;
        return res;
    }
    double *k_3(Object &object, int i){
        static double res[2] = {0, 0};
        Object object2;
        object2.x[i] = object.x[i] + 1.0/2.0 * k_2(object, i)[0];
        object2.v[i] = object.v[i] + 1.0/2.0 * k_2(object, i)[1];
        res[0] = f(object2, i)[0]*dt;
        res[1] = f(object2, i)[1]*dt;
        return res;
    }
    double *k_4(Object &object, int i){
        static double res[2] = {0, 0};
        Object object2;
        object2.x[i] = object.x[i] + k_3(object, i)[0];
        object2.v[i] = object.v[i] + k_3(object, i)[1];
        res[0] = f(object2, i)[0]*dt;
        res[1] = f(object2, i)[1]*dt;
        return res;
    }

    void count(Object &object){  
        for (int i = 0; i < n; i++){
            object.x[i+1] = object.x[i] + 1.0/6.0 * (k_1(object, i)[0] + 2*k_2(object, i)[0] + 2*k_3(object, i)[0] + k_4(object, i)[0]); 
            object.v[i+1] = object.v[i] + 1.0/6.0 * (k_1(object, i)[1] + 2*k_2(object, i)[1] + 2*k_3(object, i)[1] + k_4(object, i)[1]);
        }
    }

    RungeKutta& operator=(RungeKutta const &src) = delete;

    RungeKutta(RungeKutta const &src) = delete;
};


//---------------------------------------//
//-------SCRIPTS FOR DATA OUTPUT---------//
//---------------------------------------//
class FileOutput{
public:
    void write(Object &object){
        ofstream out;
        out.open("/home/starman/CLionProjects/RangiCut/preciseSolution.txt");
        if (out.is_open())
        {
            for (int i = 0; i < n; i++){
                out <<  dt * i  << ' ' << object.x_0*cos(object.w*dt*i) << ' ' <<  object.x_0 * sin(object.w*dt*i) * object.x_0 *sin(object.w*dt*i) / 2 + object.w * object.w *
                object.x_0 * cos(object.w*dt*i) * object.x_0 * cos(object.w*dt*i) / 2 << ' ' << -sin(object.w*dt*i) <<'\n';
            }
        }
        out.close();

        out.open("/home/starman/CLionProjects/RangiCut/RungeKuttData.txt");
        if (out.is_open())
        {
            for (int i = 0; i < n; i++){
                    out <<  dt * i  << ' ' << object.x[i]<< ' ' << object.energy[i] << ' ' << object.v[i] <<'\n';
            }
        }
        out.close();
    }
};

//-------VARIOUS BASIC TESTS AND EXPERIMANTS---------//
//--Simple Method Test------------calculate trajectory, using object's start conditions
//--KahanError--------------------estimates error that occures depending on machine epsilon
//--Time Reverse------------------calculates trajectory in strait direction and then backwards
//--Time Reverse Error------------estimates error depending on scale of dt 
class Tests{
public:


    void SimpleMethodTest(){
        Object object;
        RungeKutta method;
        method.count(object);
        object.measureEnergy();
        FileOutput file;
        file.write(object);
    }


    void KahanError() {
        ofstream out;
        out.open("/home/starman/CLionProjects/RangiCut/TimeError.txt");
        double dt_loc = dt;
        int nn = 500;
        for (int i = 0; i < nn - 2; i++) {
            double id = i;
            dt_loc = dt - dt * id / nn;
            n = floor(duration / dt_loc);
            Object object1, object2, object3;
            HoynsScheme hoynsScheme;
            HoynsSchemeKahan hoynsSchemeKahan;
            hoynsSchemeKahan.count(object1);
            hoynsScheme.count(object2);
            double nd = n;
            for (int i = 0; i < n; i++) {
                object3.x[i] = abs(object1.x[i] - object2.x[i]);
                //object3.v[i] = object1.v[i] - object2.v[i];
                //object3.energy[i] = object1.energy[i] - object2.energy[i];
            }    
            if (out.is_open()) {
                out << dt_loc << ' ' << KahanSum(object3.x, n) / nd << '\n';
            }
        }
        out.close();
    }


        void TimeReverse(){
            Tests test;
            test.SimpleMethodTest();
            Object object, object2;

            HoynsSchemeKahan method;
            method.count(object);
            object.measureEnergy();
            dt = -dt;
            object2.x[0] = object.x[n-1];
            object2.v[0] = object.v[n-1];
            method.count(object2);
            object2.measureEnergy();


            ofstream out;
            out.open("/home/starman/CLionProjects/RangiCut/timeReverse.txt");
            if (out.is_open())
            {
                for (int i = 0; i < n; i++){
                    out <<  -dt * (n-i-1)  << ' ' << object2.x[i]<< ' ' << object2.energy[i] << ' ' << object2.v[i] <<'\n';
                }
            }
            out.close();
            dt = abs(dt);
        }

        
        void TimeReverseError(){
        ofstream out;
        out.open("/home/starman/CLionProjects/RangiCut/TimeReverseError.txt");
        double dt_loc = dt;
        int nn = 500;
        for (int i = 0; i < nn; i++) {
            double id = i;
            dt_loc = dt - dt * id / nn;
            n = floor(duration / dt_loc);
            Object object1, object2, object3;
            EulersMethod method;
            method.count(object1);
            object1.measureEnergy();
            dt_loc = -dt_loc;
            object2.x[0] = object1.x[n];
            object2.v[0] = object1.v[n];
            method.count(object2);
            object2.measureEnergy();
            dt_loc = -dt_loc;
            for (int j = 0; j < n; j++) {
                object3.x[j] = abs(object1.x[j] - object2.x[n-j]);
                //object3.v[i] = object1.v[i] - object2.v[i];
                object3.energy[j] = abs(object1.energy[j] - object2.energy[n-j]);
            } 
            double nd = n;   
            if (out.is_open()) {
                out << dt_loc << ' ' << KahanSum(object3.x, n) / nd <<  ' ' << KahanSum(object3.energy, n) / nd << '\n';
            }
        }
        out.close();
        n = floor(duration/dt);
        }
};

int main(int argc, char* argv[]) {
    init("config.json");
    //Tests test;
    //test.SimpleMethodTest();
    //test.TimeReverse();
    //test.KahanError();
    //test.TimeReverseError();
    //solver<PhysPendState, GenericRK<PhysPendState>>("/home/starman/CLionProjects/RangiCut/PhysPend3.txt");
    //solver<MathPendState, GenericRK<MathPendState>>("/home/starman/CLionProjects/RangiCut/MathPend3.txt");
    //solver<PhysGammaPendState, GenericRK<PhysGammaPendState>>("/home/starman/CLionProjects/RangiCut/PhysGammaPend4.txt");
    //PGP_presise_solution prsol;
    //prsol.solve("/home/starman/CLionProjects/RangiCut/GammaPres.txt");
    solver<PhysDrivenPendState, GenericRK<PhysDrivenPendState>>("/home/starman/CLionProjects/RangiCut/PhysDrivenPend3.txt");
    DO_precise_solution prsol;
    prsol.solve("/home/starman/CLionProjects/RangiCut/DO_precise.txt");
    return 0;
}