# 
# Matlab滤波器设计APP(MatlabFilterDesignApp)
> 这是一个Matlab版本滤波器设计APP，目的是为了进行无源滤波器设计，更多可以关注微信公众号"[电路设计小工具](https://mp.weixin.qq.com/s/fxfEnir-hU0YvF9_CWyI6g)".

### 如何使用(how to use)
下载整个文件包，双击eToolsFilterDesignApp.mlapp运行即可。

### 运行环境(Env)
win10, Matlab2021a


### 程序结构

更多技术细节可以参考$^{[1][2][3][4]}$

### 特性
基于Matlab的appdesign工具开发了一套Butterworth滤波器设计软件，主要特点是：
- 支持4种不同滤波器通带类型(LPF,HPF,BPF,BRF)设计
- T型和PI型结构滤波器随意切换
- 可以设置阻带衰减决定滤波器阶数
- 可以设置通带衰减来综合滤波器
- 可以随意配置负载和终端阻抗，并支持一端接载(源端电阻短路，源端电流源，终端开路，终端短路)设计
- 可以幅频响应分析、零极点分析、瞬态分析

### 运行效果(demo)

#### LPF设计举例
设计一款-3dB截止频率为1GHz，7阶低通Butterworth滤波器，输入输出阻抗为50欧姆，设计过程如下：
![LPF滤波器设计](LPF_design_1GHz.gif)
最终设计参数如下：
![LPF滤波器设计最终](LPF_7th_Design_Final.png)
零极点仿真结果：
![LPF滤波器设计最终](LPF_7th_Design_Final_PZ.png)

#### BPF设计举例
设计7阶带通Butterworth滤波器，输入阻抗为50欧姆，输出开路，并进行瞬态仿真，设计过程如下：
![BPF滤波器设计](BPF_design_1GHz.gif)
最终设计参数如下：
![BPF滤波器设计最终](BPF_7th_Design_Final.png)
瞬态仿真结果：
![BPF滤波器设计最终](BPF_7th_Design_Final_TRAN.png)
上述设计的带通滤波器只有一种综合结果，其中Tee和PI选择并没有任何作用。

程序的Matalb源码已经上传GitHub中$^{[6]}$(https://github.com/etools361/MatlabFilterDesignApp)，有兴趣的同学可以下载试用体验，当然也欢迎技术交流。


### 参考(ref)
[1] : [模拟无源滤波器设计（一）](https://mp.weixin.qq.com/s/wNRHyBHpimjU90bymHp7JA) \
[2] : [模拟无源滤波器设计（二）](https://mp.weixin.qq.com/s/3GMQs4WDm683tdAXqyoOgQ) \
[3] : [模拟无源滤波器设计（三）](https://mp.weixin.qq.com/s/nZFx7weLcO-WRKLbP0T4jQ) \
[4] : [模拟无源滤波器设计（四）](https://mp.weixin.qq.com/s/mllwGShvbh3TWdFRbp9LhQ) \
[5] : [欢迎关注电路设计小工具公众号](https://mp.weixin.qq.com/s/fxfEnir-hU0YvF9_CWyI6g) \
[6] : [Matlab源码地址](https://github.com/etools361/MatlabFilterDesignApp)
