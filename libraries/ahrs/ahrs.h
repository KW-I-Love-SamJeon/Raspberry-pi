#include <math.h>

typedef struct{
    double w, x, y, z;
}Quaternion;

void QuaternInsert(Quaternion *a, double w, double x, double y, double z){
    a->w = w;
    a->x = x;
    a->y = y;
    a->z = z;
}

void QuaternProd(Quaternion *c, Quaternion a, Quaternion b){
    c->w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
    c->x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y;
    c->y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x;
    c->z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w;
};

void QuaternScalar(Quaternion *b, Quaternion a, double s){
    b->w = a.w * s;
    b->x = a.x * s;
    b->y = a.y * s;
    b->z = a.z * s;
}

void QuaternAdd(Quaternion *c, Quaternion a, Quaternion b){
    c->w = a.w + b.w;
    c->x = a.x + b.x;
    c->y = a.y + b.y;
    c->z = a.z + b.z; 
}

void QuaternSub(Quaternion *c, Quaternion a, Quaternion b){
    c->w = a.w - b.w;
    c->x = a.x - b.x;
    c->y = a.y - b.y;
    c->z = a.z - b.z; 
}

Quaternion QuaternConj(Quaternion a){
    a.x = -a.x;
    a.y = -a.y;
    a.z = -a.z;
    return a;
}

// SE_q: quaternion 1X4
// E_d: 1X4 reference 물리량
// S_s: 1X4 측정 물리량
// output : 3X1
void GradientDescentLoss(double (*output)[1], Quaternion SE_q
                        , Quaternion E_d, Quaternion S_s){

    output[0][0] = 2.0 * E_d.x * (0.5 - pow(SE_q.y, 2) - pow(SE_q.z, 2)) 
                + 2.0 * E_d.y * (SE_q.w * SE_q.z + SE_q.x * SE_q.y)
                + 2.0 * E_d.z * (SE_q.x * SE_q.z - SE_q.w * SE_q.y) - S_s.x;

    output[1][0] = 2.0 * E_d.x * (SE_q.x * SE_q.y - SE_q.w * SE_q.z)
                + 2.0 * E_d.y * (0.5 - pow(SE_q.x, 2) - pow(SE_q.z, 2))
                + 2.0 * E_d.z * (SE_q.w * SE_q.x + SE_q.y * SE_q.z) - S_s.y;

    output[2][0] = 2.0 * E_d.x * (SE_q.w * SE_q.y + SE_q.x * SE_q.z)
                + 2.0 * E_d.y * (SE_q.y * SE_q.z - SE_q.w * SE_q.x)
                + 2.0 * E_d.z * (0.5 - pow(SE_q.x, 2) - pow(SE_q.y, 2)) - S_s.z;
}

// output : quaternion 3X4 
void JacobianLoss(double (*output)[4], Quaternion SE_q, Quaternion E_d){

    output[0][0] = 2.0 * (E_d.y * SE_q.z - E_d.z * SE_q.y);
    output[0][1] = 2.0 * (E_d.y * SE_q.y - E_d.z * SE_q.z);
    output[0][2] = -4.0 * E_d.x * SE_q.y + 2.0 * (E_d.y * SE_q.x - E_d.z * SE_q.w);
    output[0][3] = -4.0 * E_d.x * SE_q.z + 2.0 * (E_d.y * SE_q.w + E_d.z * SE_q.x);

    output[1][0] = 2.0 * (-E_d.x *SE_q.z + E_d.z * SE_q.x);
    output[1][1] = -4.0 * E_d.y * SE_q.x + 2.0 * (E_d.x * SE_q.y + E_d.z * SE_q.w);
    output[1][2] = 2.0 * (E_d.x * SE_q.x + E_d.z * SE_q.z);
    output[1][3] = -4.0 * E_d.y * SE_q.z + 2.0 * (-E_d.x * SE_q.w + E_d.z * SE_q.y);

    output[2][0] = 2.0 * (E_d.x * SE_q.y - E_d.y * E_d.x);
    output[2][1] = -4.0 * E_d.z * SE_q.x + 2.0 * (E_d.x * SE_q.z - E_d.y * SE_q.w);
    output[2][2] = -4.0 * E_d.z * SE_q.y + 2.0 * (E_d.x * SE_q.w + E_d.y * SE_q.z);
    output[2][3] = 2.0 * (E_d.x * SE_q.x + E_d.y * E_d.y);
}

void StepCal(Quaternion *output, double (*f)[1], double (*d)[4]){
    double tp_d[4][6];
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 6; j++)
            tp_d[i][j] = d[j][i];
    
    double matmul[4][1] = { 0 };
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 6; j++) 
            matmul[i][0] += tp_d[i][j] * f[j][0];
    
    double step_norm = sqrt(pow(matmul[0][0], 2) + pow(matmul[1][0], 2)
                     + pow(matmul[2][0], 2) + pow(matmul[3][0], 2));
    
    QuaternInsert(output, matmul[0][0], matmul[1][0], matmul[2][0], matmul[3][0]);
    QuaternScalar(output, *output, 1 / step_norm);
}

void Madgwick_9_DOF(Quaternion *output, Quaternion q, Quaternion Accelerometer
                    , Quaternion Gyroscope, Quaternion Magnetometer, double beta, double period){

    double f[6][1], d[6][4];

    // Accelerometer Measurement
    double accel_norm = sqrt(pow(Accelerometer.x, 2) + pow(Accelerometer.y, 2) + pow(Accelerometer.z, 2));
    Quaternion accel_normalized;
    QuaternScalar(&accel_normalized, Accelerometer, 1 / accel_norm);

    Quaternion r_vector;
    QuaternInsert(&r_vector, 0, 0, 0, 1);
    GradientDescentLoss(f, q, r_vector, accel_normalized);
    JacobianLoss(d, q, r_vector);


    // Magnetometer Measurement
    double magnet_norm = sqrt(pow(Magnetometer.x, 2) + pow(Magnetometer.y, 2) + pow(Magnetometer.z, 2));
    Quaternion magnet_normalized;
    QuaternScalar(&magnet_normalized, Magnetometer, 1 / magnet_norm); 

    Quaternion e_magnet;
    QuaternProd(&e_magnet, magnet_normalized, QuaternConj(q));
    QuaternProd(&e_magnet, q, e_magnet);
    QuaternInsert(&e_magnet, 0, sqrt(pow(e_magnet.x, 2) + pow(e_magnet.y, 2)), 0, e_magnet.z);

    GradientDescentLoss(f + 3, q, e_magnet, magnet_normalized);
    JacobianLoss(d + 3, q, e_magnet);

    Quaternion step;
    StepCal(&step, f, d);
    
    Quaternion q_dot;
    QuaternProd(&q_dot, q, Gyroscope);
    QuaternScalar(&q_dot, q_dot, 0.5);
    QuaternScalar(&step, step, beta);
    QuaternSub(&q_dot, q_dot, step);

    QuaternScalar(&q_dot, q_dot, period);
    QuaternAdd(output, q, q_dot);
    double q_norm = sqrt(pow(output->w, 2) + pow(output->x, 2)
                  + pow(output->y, 2) + pow(output->z, 2));
    QuaternScalar(output, *output, 1 / q_norm);
}