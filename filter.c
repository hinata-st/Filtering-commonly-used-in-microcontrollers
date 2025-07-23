/*
 * filter.c
 *
 *  Created on: 2022年12月27日
 *      Author: Emilia-tan
 */

/***************************************************
函数名: float kalman_filter(float raw)
说明:	非矩阵一维卡尔曼滤波
入口:	float raw-需要滤波的值 
出口:	x_t - 滤波后的值 
备注:	Q,R 属于超参数，需要手动调节
		k = ((p_t + Q)/(p_t + Q + R)) k大相信估计值，k小相信测量值 
		此滤波适用于  状态量 = 观测量 
****************************************************/
#include "filter.h"
#include "math.h"


int text(one_struct aaa)
{
	aaa.a = 1;
}


typedef struct {
    float X_last; //上一时刻的最优结果  X(k-|k-1)
    float X_mid;  //当前时刻的预测结果  X(k|k-1)
    float X_now;  //当前时刻的最优结果  X(k|k)
    float P_mid;  //当前时刻预测结果的协方差  P(k|k-1)
    float P_now;  //当前时刻最优结果的协方差  P(k|k)
    float P_last; //上一时刻最优结果的协方差  P(k-1|k-1)
    float kg;     //kalman增益
    float A;      //系统参数
    float B;
    float Q;
    float R;
    float H;
}extKalman_t;


/**
  * @name   kalmanCreate
  * @brief  创建一个卡尔曼滤波器
  * @param  p:  滤波器
  *         T_Q:系统噪声协方差
  *         T_R:测量噪声协方差
  *         
  * @retval none
  * @attention R固定，Q越大，代表越信任侧量值，Q无穷代表只用测量值
  *		       	反之，Q越小代表越信任模型预测值，Q为零则是只用模型预测
  */
void KalmanCreate(extKalman_t *p,float T_Q,float T_R)
{
    p->X_last = (float)0;
    p->P_last = 0;
    p->Q = T_Q;
    p->R = T_R;
    p->A = 1;
    p->B = 0;
    p->H = 1;
    p->X_mid = p->X_last;
}

/**
  * @name   KalmanFilter
  * @brief  卡尔曼滤波器
  * @param  p:  滤波器
  *         dat:待滤波数据
  * @retval 滤波后的数据
  * @attention Z(k)是系统输入,即测量值   X(k|k)是卡尔曼滤波后的值,即最终输出
  *            A=1 B=0 H=1 I=1  W(K)  V(k)是高斯白噪声,叠加在测量值上了,可以不用管
  *            以下是卡尔曼的5个核心公式
  *            一阶H'即为它本身,否则为转置矩阵
  */

float KalmanFilter(extKalman_t* p,float dat)
{
    p->X_mid =p->A*p->X_last;                     //百度对应公式(1)    x(k|k-1) = A*X(k-1|k-1)+B*U(k)+W(K)     状态方程
    p->P_mid = p->A*p->P_last+p->Q;               //百度对应公式(2)    p(k|k-1) = A*p(k-1|k-1)*A'+Q            观测方程
    p->kg = p->P_mid/(p->P_mid+p->R);             //百度对应公式(4)    kg(k) = p(k|k-1)*H'/(H*p(k|k-1)*H'+R)   更新卡尔曼增益
    p->X_now = p->X_mid + p->kg*(dat-p->X_mid);   //百度对应公式(3)    x(k|k) = X(k|k-1)+kg(k)*(Z(k)-H*X(k|k-1))  修正估计值
    p->P_now = (1-p->kg)*p->P_mid;                //百度对应公式(5)    p(k|k) = (I-kg(k)*H)*P(k|k-1)           更新后验估计协方差
    p->P_last = p->P_now;                         //状态更新
    p->X_last = p->X_now;
    return p->X_now;							  //输出预测结果x(k|k)
}


float kalman_filter_liner(float raw)
{
	static char init = 0;
	//估计值 
	static float x_t = 0;
	//估计协方差 
	static float p_t = 0;
	//状态转移协方差 
	static float Q = 1;
	//测量噪声协方差 
	static float R = 30;
	//卡尔曼增益 
	static float k = 0;
	
	//是否第一次运行 
	if(init == 0)
	{
		//进行参数初始化
		x_t = 0;
		p_t = 1;
	}
	//预测 
	//状态转移方程x(t) = x(t-1) 
	//x_t = x_t;
	p_t = p_t + Q;
	//更新 
	k = p_t/(p_t+R);
	x_t = x_t + k*(raw - x_t);
	p_t = (1 - k)*p_t;	
	return x_t;
}


//核心思想：根据当前的仪器"测量值" 和上一刻的 “预测量” 和 “误差”，计算得到当前的最优量，再预测下一刻的量。里面比较突出的是观点是：把误差纳入计算，而且分为预测误差和测量误差两种，通称为噪声。
//还有一个非常大的特点是：误差独立存在，始终不受测量数据的影响。
//优点：巧妙的融合了观测数据与估计数据，对误差进行闭环管理，将误差限定在一定范围。适用性范围很广，时效性和效果都很优秀。
//缺点：需要调参，参数的大小对滤波的效果影响较大。
//卡尔曼滤波
int KalmanFilter(int inData)
{
      static float prevData = 0;                                 //先前数值
      static float p = 10, q = 0.001, r = 0.001, kGain = 0;      // q控制误差  r控制响应速度 

      p = p + q;
      kGain = p / ( p + r );                                     //计算卡尔曼增益
      inData = prevData + ( kGain * ( inData - prevData ) );     //计算本次滤波估计值
      p = ( 1 - kGain ) * p;                                     //更新测量方差
      prevData = inData;
      return inData;                                             //返回滤波值
}

//-------------------------------------------------------------------------------------------------------------------
// 函数简介     卡尔曼滤波器陀螺仪滤波 
// 参数说明     float angle_m(陀螺仪加速度计估计角度),float gyro_m(陀螺仪角速度数据)
// 返回参数     angle_l (滤波后的角度)
// 使用示例     Kalman_Pose_Y((atan2f(acc_x,acc_z)*180/3.1416),gyro_y,0.001);
//				Kalman_Pose_X((atan2f(acc_y,acc_z)*180/3.1416), gyro_x,0.001);
//				Q_angle,Q_gyro,R_angle属于超参数，需要手动调节 
// 备注信息     无
//-------------------------------------------------------------------------------------------------------------------
float kalman_filter(float angle_m,float gyro_m,float dt)          //gyro_m:gyro_measure
{
	//Q_angle(角度噪声)  Q_gyro(角速度漂移噪声) R_angle(角度测量噪声)
    static float Q_angle=0.001, Q_gyro=0.003, R_angle=0.03,angle_l = 0;       //注意：dt的取值为kalman滤波器采样时间;
    static float P[2][2] = {
                                { 1, 0 },
                                { 0, 1 }
                            };
    static float Pdot[4] ={0,0,0,0};
    static const char C_0 = 1;
    static float q_bias, angle_err, PCt_0, PCt_1, E, kk_0, kk_1, t_0, t_1,flag = 0;
    //如果是第一次运行 
    if(flag == 0)
    {
		flag = 1;
		angle_l = angle_m;
	}
    //1时间更新(预测): X(k|k-1) = A(k,k-1)*X(k-1|k-1) + B(k)*u(k)
    angle_l = angle_l + (gyro_m-q_bias) * dt;
	//2更新先验协方差: P(k|k-1) = A(k,k-1)*A(k,k-1)^T*P(k-1|k-1)+Q(k)
    Pdot[0]= - P[0][1] - P[1][0];
    Pdot[1]=- P[1][1];
    Pdot[2]=- P[1][1];
    Pdot[3]=Q_gyro;

    P[0][0] += Pdot[0] * dt	+ Q_angle;
    P[0][1] += Pdot[1] * dt;
    P[1][0] += Pdot[2] * dt;
    P[1][1] += Pdot[3];
    //3测量方程:Z(k) = HX(k) + V(K) 
    
	//4计算卡尔曼增益: K(k) = P(k|k-1)*H(k)^T/(P(k|k-1)*H(k)*H(k)^T + R(k))
    PCt_0 = C_0 * P[0][0];
    PCt_1 = C_0 * P[1][0];
    E = R_angle + C_0 * PCt_0;
    kk_0 = PCt_0 / E;
    kk_1 = PCt_1 / E;    
	//5测量更新(校正): X(k|k) = X(k|k-1)+K(k)*(Z(k)-H(k)*X(k|k-1))
	angle_err = angle_m - angle_l;	
    angle_l += kk_0 * angle_err;
    q_bias  += kk_1 * angle_err;          
	//6更新后验协方差: P(k|k) =（I-K(k)*H(k))*P(k|k-1)
    t_0 = PCt_0;
    t_1 = C_0 * P[0][1];
    P[0][0] -= kk_0 * t_0;
    P[0][1] -= kk_0 * t_1;
    P[1][0] -= kk_1 * t_0;
    P[1][1] -= kk_1 * t_1;

    return angle_l;
}


/***************************************************
函数名: float rc_filter_lon(float data_now,float p)
说明:	一阶RC低通滤波 
入口:	float data_now-需要滤波的值 
出口:	data_filter - 滤波后的值 
备注:	p 属于超参数，需要手动调节 
		p大相信当前值，p小相信过去值 
****************************************************/
float rc_filter_lon(float data_now,float p)
{
    static float last_data;
    static float data_filter;
    static char flag;

    if(flag == 0)//如果是第一次运行
    {
        last_data = data_now;
        flag++;
    }

    data_filter = data_now*p + last_data*(1-p);

    last_data = data_filter;

    return data_filter;
}


/***************************************************
函数名: float filter(float current)
说明:	滑动滤波 
入口:	float current-需要滤波的值 
出口:	sum/K_1 - 滤波后的值 
		0xFFFF - 当前值无效 
备注:	SUM_WIN_SIZE滑动窗口大小 
****************************************************/
#define SUM_WIN_SIZE  6
float filter(float current)
{
	static float history[SUM_WIN_SIZE];//历史值，其中history[SUM_WIN_SIZE-1]为最近的记录
	static int buff_init=0;//前SUM_WIN_SIZE-1次填充后才能开始输出
	static int index1=0;//环形数组可放数据的位置
	static int factor[SUM_WIN_SIZE]={1,2,3,4,5,6};//加权系数
	int K_1=21;//1+2+3+4+5+6	
		
    int i,j;
    float sum=0;

    if(buff_init==0)
    {
        history[index1]=current;
        index1++;
        if(index1>=(SUM_WIN_SIZE-1))
        {
            buff_init=1;//index有效范围是0-5，前面放到5，下一个就可以输出
        }
        return 0xFFFF;//当前无法输出，做个特殊标记区分
    }
    else
    {
        history[index1]=current;
        index1++;
        if(index1>=SUM_WIN_SIZE)
        {
            index1=0;//index有效最大5,下次再从0开始循环覆盖
        }

        j=index1;
        for(i=0;i<SUM_WIN_SIZE;i++)
        {
            //注意i=0的值并不是最早的值
            sum+=history[j]*factor[i];//注意防止数据溢出
            j++;
            if(j==SUM_WIN_SIZE)
            {
                j=0;
            }
        }
        return sum/K_1;
    }
}

/***************************************************
函数名: float LPButterworth(float curr_input,lpf_buf *buf,lpf_param *params)
说明:	二阶巴特沃斯低通滤波器
入口:	float curr_input-当前滤波器输入
			lpf_buf *buf-滤波器中间状态
			lpf_param *params-滤波器参数
出口:	float 滤波器输出值
备注:	无
****************************************************/
typedef struct
{
 float input[3];
 float output[3];
}lpf_buf;

typedef struct
{
  float a[3];
  float b[3];
}lpf_param;

float LPButterworth(float curr_input,lpf_buf *buf,lpf_param *params)
{
	if(buf->output[0]==0&&
		 buf->output[1]==0&&
		 buf->output[2]==0&&
		 buf->input[0]==0&&
		 buf->input[1]==0&&
		 buf->input[2]==0)
	{
		buf->output[0]=curr_input;
		buf->output[1]=curr_input;
		buf->output[2]=curr_input;
		buf->input[0]=curr_input;
		buf->input[1]=curr_input;
		buf->input[2]=curr_input;
		return curr_input;
	}
	
  /* 加速度计Butterworth滤波 */
  /* 获取最新x(n) */
  buf->input[2]=curr_input;
  /* Butterworth滤波 */
  buf->output[2]=params->b[0] * buf->input[2]
													+params->b[1] * buf->input[1]
													+params->b[2] * buf->input[0]
													-params->a[1] * buf->output[1]
													-params->a[2] * buf->output[0];
  /* x(n) 序列保存 */
  buf->input[0]=buf->input[1];
  buf->input[1]=buf->input[2];
  /* y(n) 序列保存 */
  buf->output[0]=buf->output[1];
  buf->output[1]=buf->output[2];
	
	for(int i=0;i<3;i++)
	{
	  if(isnan(buf->output[i])==1
			||isnan(buf->input[i])==1)		
			{		
				buf->output[0]=curr_input;
				buf->output[1]=curr_input;
				buf->output[2]=curr_input;
				buf->input[0]=curr_input;
				buf->input[1]=curr_input;
				buf->input[2]=curr_input;
				return curr_input;
			}
	}	
  return buf->output[2];
}

/***************************************************
函数名: void set_cutoff_frequency(float sample_frequent, float cutoff_frequent,lpf_param *LPF)
说明:	二阶巴特沃斯低通滤波器参数设计
入口:	float sample_frequent-采样频率
			float cutoff_frequent-截止频率
			lpf_param *LPF-滤波器参数
出口:	无
备注:	无
****************************************************/
#define M_PI_F 3.141592653589793f
void set_cutoff_frequency(float sample_frequent, float cutoff_frequent,lpf_param *LPF)
{
  float fr = sample_frequent / cutoff_frequent;
  float ohm = tanf(M_PI_F / fr);
  float c = 1.0f + 2.0f * cosf(M_PI_F / 4.0f) * ohm + ohm * ohm;
  if (cutoff_frequent <= 0.0f) {
    // no filtering
    return;
  }
  LPF->b[0] = ohm * ohm / c;
  LPF->b[1] = 2.0f * LPF->b[0];
  LPF->b[2] = LPF->b[0];
  LPF->a[0]=1.0f;
  LPF->a[1] = 2.0f * (ohm * ohm - 1.0f) / c;
  LPF->a[2] = (1.0f - 2.0f * cosf(M_PI_F / 4.0f) * ohm + ohm * ohm) / c;
}

//用法示例 
//lpf_param accel_lpf_param,gyro_lpf_param;
//lpf_buf gyro_filter_buf[3],accel_filter_buf[3];
//set_cutoff_frequency(200,20,&accel_lpf_param);
//smartcar_imu.accel_g_raw.z=LPButterworth(smartcar_imu._accel_g_raw.z,&accel_filter_buf[2],&accel_lpf_param);


/***************************************************
函数名: float LowOassFilter(float raw)
说明:	自适应低通滤波 
入口:	float raw-需要滤波的数据 
		float Tf- 低通滤波时间常数 
出口:	y - 滤波器输出值 
备注:	无
****************************************************/
float LowOassFilter(float raw,float Tf)
{
	//最后执行时间戳 
	static unsigned long  timestamp_prev = 0;
	//当前执行时间戳
	static unsigned long timestamp = 0;
	//上一个循环中的过滤后的值 
	static unsigned long y_prev = 0;
	//时间间隔 
	static float dt = 0;
	//是否第一次执行 
	static char flag = 0; 
	//需要自己实现获取时间戳 
    //unsigned long timestamp = micros();
    if(flag == 0)//第一次执行 
    {
		timestamp_prev = timestamp;
		y_prev =  raw;
		flag = 1;
		return  y_prev;
	}
	else
	{
		dt = (timestamp - timestamp_prev)*1e-6f;
	}
    if (dt < 0.0f ) dt = 1e-3f;
    else if(dt > 0.3f) {
        y_prev = raw;
        timestamp_prev = timestamp;
        return raw;
    }
    float alpha = Tf/(Tf + dt);
    float y = alpha*y_prev + (1.0f - alpha)*raw;
    y_prev = y;
    timestamp_prev = timestamp;
    return y;
}

//方法：连续采样N次（N取奇数）把N次采样值按大小排列取中间值为本次有效值
//优点：能有效克服因偶然因素引起的波动干扰；对温度、液位等变化缓慢的被测参数有良好的滤波效果
//缺点：对流量，速度等快速变化的参数不宜
//中值滤波算法
int middleValueFilter(int N,int data)
{
      int value_buf[N];
      int i,j,k,temp;
      for( i = 0; i < N; ++i)
      {
        value_buf[i] = data;  

      }

      for(j = 0 ; j < N-1; ++j)
      {
          for(k = 0; k < N-j-1; ++k)
          {
              //从小到大排序，冒泡法排序
              if(value_buf[k] > value_buf[k+1])
              {
                temp = value_buf[k];
                value_buf[k] = value_buf[k+1];
                value_buf[k+1] = temp;
              }
          }
      }

      return value_buf[(N-1)/2];
}


//方法：连续取N个采样值进行算术平均运算;
//N值较大时：信号平滑度较高，但灵敏度较低
//N值较小时：信号平滑度较低，但灵敏度较高
//N值的选取：一般流量，N=12；压力：N=4
//优点：试用于对一般具有随机干扰的信号进行滤波。这种信号的特点是有一个平均值，信号在某一数值范围附近上下波动。
//缺点：测量速度较慢或要求数据计算较快的实时控制不适用。
//算术平均值滤波
int averageFilter(int N,int data)
{
     int sum = 0;
     short i;
     for(i = 0; i < N; ++i)
     {
        sum += data;  
     }
     return sum/N;
}



//方法：把连续取N个采样值看成一个队列，队列的长度固定为N。每次采样到一个新数据放入队尾，并扔掉原来队首的一次数据(先进先出原则)。把队列中的N个数据进行算术平均运算,就可获得新的滤波结果。
//N值的选取：流量，N=12；压力：N=4；液面，N=4~12；温度，N=1~4
//优点：对周期性干扰有良好的抑制作用，平滑度高；试用于高频振荡的系统
//缺点：灵敏度低；对偶然出现的脉冲性干扰的抑制作用较差，不适于脉冲干扰较严重的场合
//比较浪费RAM（改进方法，减去的不是队首的值，而是上一次得到的平均值）
//平滑均值滤波
#define N 10
int value_buf[N];
int sum=0;
int curNum=0;

int moveAverageFilter(int data)
{
      if(curNum < N)
      {
          value_buf[curNum] = data;
          sum += value_buf[curNum];
          curNum++;
          return sum/curNum;
      }
      else
      {
          sum -= sum/N;
          sum += data;
          return sum/N;
      }
}


//方法：相当于“限幅滤波法”+“递推平均滤波法”
//每次采样到的新数据先进行限幅处理再送入队列进行递推平均滤波处理
//优点：对于偶然出现的脉冲性干扰，可消除有其引起的采样值偏差。
//缺点：比较浪费RAM
//限幅平均滤波
#define A 50        //限制幅度阈值
#define M 12
int data[M];
int First_flag=0;

int LAverageFilter(int data_new)
{
    int i;
    int temp,sum,flag=0;
    data[0] = data_new;
    for(i=1;i<M;i++)
    {
      temp=data_new;
      if((temp-data[i-1])>A || ((data[i-1]-temp)>A))
      {
          i--;
          flag++;
      }
      else
      {
          data[i]=temp;
      }
    }

    for(i=0;i<M;i++)
    {
      sum+=data[i];
    } 
    return  sum/M;
}


int main()
{
	return 0;
}
