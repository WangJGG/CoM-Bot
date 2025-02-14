# CoM-Bot

## Authors:  
lyn, wjg, zyn, wsy, byl

## Description  
CoM-Bot is a software system designed to track and train human Center of Mass (CoM) movements using advanced techniques in motion tracking and machine learning. The system integrates various tools to achieve accurate real-time tracking and prediction of a person's balance and posture. CoM-Bot utilizes MediaPipe for human pose detection, OpenSim for biomechanical modeling, LSTM networks for center of mass prediction, and Flask for building a simple web interface.

## Tech Stack  
- **MediaPipe**: A framework for building multimodal applied machine learning pipelines, used here for human pose and keypoint detection.
- **OpenSim**: A biomechanical modeling and simulation software, employed for modeling the human body and generating joint angles for center of mass estimation.
- **LSTM (Long Short-Term Memory)**: A type of recurrent neural network used for predicting the center of mass based on input features such as foot pressure sensors.
- **Flask**: A lightweight web framework for serving the project’s functionalities through a simple and interactive user interface.

## Project Structure  
- **demo.py & template**: A simple demonstration of the system’s core functionalities. This serves as an initial showcase for how the system works.
- **PiezoElectric**: Implements a system that uses foot pressure sensors to detect pressure changes and predicts the user's center of mass using an LSTM model.
- **ViconMarkerless2Opensim**: Utilizes MediaPipe to detect human body keypoints, which are then converted into joint angles in the OpenSim model to track the user’s center of mass.
