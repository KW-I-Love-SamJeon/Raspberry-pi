#include <ahrs.h>
#include <Waveshare_10Dof-D.h>
#include "Wire.h"
#include "I2Cdev.h"
#include "MPU9250.h"

MPU9250 accelgyro;
I2Cdev I2C_M;
uint8_t buffer_m[6];

int16_t ax, ay, az;
int16_t gx, gy, gz;
int16_t mx, my, mz;
int16_t qw, qx, qy, qz;

Quaternion Accelerometer;
Quaternion Gyroscope;
Quaternion Magnetometer;
Quaternion q;

float heading;
float tiltheading;

float pitch;
float roll;

//Angle
float a_pitch;
float a_roll;

#define sample_num_mdate 5000
volatile float mx_sample[3];  
volatile float my_sample[3];
volatile float mz_sample[3];
static float mx_centre = 0;
static float my_centre = 0;
static float mz_centre = 0;
volatile int mx_max = 0;
volatile int my_max = 0;
volatile int mz_max = 0;
volatile int mx_min = 0;
volatile int my_min = 0;
volatile int mz_min = 0;

unsigned long readTime; // timer

unsigned long time_previous; // period
unsigned long time_current;
double period;

void setup()
{
Wire.begin();
Serial.begin(9600);

Serial.println("CLEARDATA"); // clear excel data
Serial.println("LABEL, #timestamp [ns], p_RS_R_x [m], p_RS_R_y [m], p_RS_R_z [m], q_RS_w [], q_RS_x [], q_RS_y [], q_RS_z [], v_RS_R_x [m s^-1], v_RS_R_y [m s^-1], v_RS_R_z [m s^-1],b_w_RS_S_x [rad s^-1], b_w_RS_S_y [rad s^-1], b_w_RS_S_z [rad s^-1], b_a_RS_S_x [m s^-2], b_a_RS_S_y [m s^-2], b_a_RS_S_z [m s^-2], pitch, roll, yaw");

accelgyro.initialize();
delay(1000);
time_previous = millis();
}
void loop()
{
get_Data();
getCompassDate_calibrated(); 
getHeading(); 
getTiltHeading();
getAngle();

//Time
Serial.print("DATA,");

readTime = millis(); // ms
readTime = readTime * 1000 * 1000; // ns 
Serial.print(readTime);
Serial.print(",");


// Magnetic
Serial.print(Magnetometer.x,16);
Serial.print(",");
Serial.print(Magnetometer.y,16);
Serial.print(",");
Serial.print(Magnetometer.z,16);
Serial.print(",");
//

//Quaternion
QuaternInsert(&q, 1, 0, 0, 0);
time_current = millis(); // ms
period = double(time_current - time_previous) / 1000; // sec
time_previous = time_current;

Madgwick_9_DOF(&q, q,
               Accelerometer,
               Gyroscope,
               Magnetometer,
               0.05, period);

Serial.print(q.w, 16);
Serial.print(",");
Serial.print(q.x, 16);
Serial.print(",");
Serial.print(q.y, 16);
Serial.print(",");
Serial.print(q.z, 16);
Serial.print(",");
//

// Velocity (test zero)
Serial.print(mx_centre);
Serial.print(",");
Serial.print(my_centre);
Serial.print(",");
Serial.print(mz_centre);
Serial.print(",");
//

// Gyro
Serial.print(Gyroscope.x,16);
Serial.print(",");
Serial.print(Gyroscope.y,16);
Serial.print(",");
Serial.print(Gyroscope.z,16);
Serial.print(",");
//

// Acc
Serial.print(Accelerometer.x,16);
Serial.print(",");
Serial.print(Accelerometer.y,16);
Serial.print(",");
Serial.print(Accelerometer.z,16);
Serial.print(",");

//
Serial.print(a_pitch);
Serial.print(",");
Serial.print(a_roll);
Serial.print(",");
Serial.print(heading); // yaw
Serial.print(",");

//
Serial.println(" ");
Serial.println();

delay(5); // 0.005sec
}


void getAngle(void)
{
  a_pitch = asin(-Accelerometer.x);
  a_roll = asin(Accelerometer.y / cos(a_pitch)) * 180 / PI;
  a_pitch = a_pitch * 180 / PI;
  if (a_pitch < 0) a_pitch += 360;
  if (a_roll < 0) a_roll += 360;
}

void getHeading(void)
{
heading = 180 * atan2(Magnetometer.y, Magnetometer.x) / PI;
if (heading < 0) heading += 360;
}

void getTiltHeading(void)
{
pitch = asin(-Accelerometer.x);
roll = asin(Accelerometer.y / cos(pitch));
float xh = Magnetometer.x * cos(pitch) + Magnetometer.z * sin(pitch);
float yh = Magnetometer.x * sin(roll) * sin(pitch) + Magnetometer.y * cos(roll) - Magnetometer.z * sin(roll) * cos(pitch);
float zh = -Magnetometer.x * cos(roll) * sin(pitch) + Magnetometer.y * sin(roll) + Magnetometer.z * cos(roll) * cos(pitch);
tiltheading = 180 * atan2(yh, xh) / PI;
if (yh < 0) tiltheading += 360;
}

void Mxyz_init_calibrated ()
{
while (!Serial.find("ready"));
get_calibration_Data ();
}

void get_calibration_Data ()
{
for (int i = 0; i < sample_num_mdate; i++)
{
get_one_sample_date_mxyz();
if (mx_sample[2] >= mx_sample[1])mx_sample[1] = mx_sample[2];
if (my_sample[2] >= my_sample[1])my_sample[1] = my_sample[2]; //find max value
if (mz_sample[2] >= mz_sample[1])mz_sample[1] = mz_sample[2];
if (mx_sample[2] <= mx_sample[0])mx_sample[0] = mx_sample[2];
if (my_sample[2] <= my_sample[0])my_sample[0] = my_sample[2]; //find min value
if (mz_sample[2] <= mz_sample[0])mz_sample[0] = mz_sample[2];
}
mx_max = mx_sample[1];
my_max = my_sample[1];
mz_max = mz_sample[1];
mx_min = mx_sample[0];
my_min = my_sample[0];
mz_min = mz_sample[0];
mx_centre = (mx_max + mx_min) / 2;
my_centre = (my_max + my_min) / 2;
mz_centre = (mz_max + mz_min) / 2;
}
void get_one_sample_date_mxyz()
{
getCompass_Data();
mx_sample[2] = Magnetometer.x;
my_sample[2] = Magnetometer.y;
mz_sample[2] = Magnetometer.z;
}

void get_Data(void)
{
accelgyro.getMotion9(&ax, &ay, &az, &gx, &gy, &gz, &mx, &my, &mz);
QuaternInsert(&Accelerometer, 0, double(ax) / 16384, double(ay) / 16384, double(az) / 16384);
QuaternInsert(&Gyroscope, 0, double(gx) * 250 / 32768, double(gy) * 250 / 32768, double(gz) * 250 / 32768);
}

void getCompass_Data(void)
{
I2C_M.writeByte(MPU9150_RA_MAG_ADDRESS, 0x0A, 0x01); //enable the magnetometer
delay(10);
I2C_M.readBytes(MPU9150_RA_MAG_ADDRESS, MPU9150_RA_MAG_XOUT_L, 6, buffer_m);
mx = ((int16_t)(buffer_m[1]) << 8) | buffer_m[0] ;
my = ((int16_t)(buffer_m[3]) << 8) | buffer_m[2] ;
mz = ((int16_t)(buffer_m[5]) << 8) | buffer_m[4] ;
QuaternInsert(&Magnetometer, 0, double(mx) * 1200 / 4096, double(my) * 1200 / 4096, double(mz) * 1200 / 4096);
}
void getCompassDate_calibrated ()
{
getCompass_Data();
Magnetometer.x = Magnetometer.x - mx_centre;
Magnetometer.y = Magnetometer.y - my_centre;
Magnetometer.z = Magnetometer.z - mz_centre;
}