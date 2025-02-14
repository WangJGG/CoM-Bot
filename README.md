# CoM-Bot
## Author: 
lyn,wjg,zyn,wsy,byl
## Description
一个人体重心跟踪与训练的软件机器人
## 技术栈
- MediaPipe
- Opensim
- LSTM
- Flask
## 项目结构
- demo.py & template: 一个简单的演示demo
- PiezoElectric: 足底压力传感，利用LSTM模型进行重心预测
- ViconMarkerless2Opensim:利用MediaPipe进行人体关键点检测，将关键点坐标转换为Opensim模型的关节角度,实现人体重心跟踪