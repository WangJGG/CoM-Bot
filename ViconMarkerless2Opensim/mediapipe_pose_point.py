#mediapipe_3D
import cv2
import mediapipe as mp
import numpy as np
import csv

mp_pose = mp.solutions.pose
mp_drawing = mp.solutions.drawing_utils
mp_drawing_styles = mp.solutions.drawing_styles
mp_holistic = mp.solutions.holistic

def extract_pose_landmarks(image, holistic):
    results = holistic.process(image)
    if results.pose_world_landmarks:
        frame_data = [[landmark.x, landmark.z, landmark.y*-1] for landmark in results.pose_world_landmarks.landmark]
        return frame_data
    return []

def draw_landmarks(image, results):
    mp_drawing.draw_landmarks(image, results.face_landmarks, mp_holistic.FACEMESH_CONTOURS, None, mp_drawing_styles.get_default_face_mesh_contours_style())
    mp_drawing.draw_landmarks(image, results.pose_landmarks, mp_holistic.POSE_CONNECTIONS, mp_drawing_styles.get_default_pose_landmarks_style())
    mp_drawing.draw_landmarks(image, results.left_hand_landmarks, mp_holistic.HAND_CONNECTIONS)
    mp_drawing.draw_landmarks(image, results.right_hand_landmarks, mp_holistic.HAND_CONNECTIONS)
    return image

def save_landmarks_to_csv(data, filename):
    all_frames_data_np = np.array(data)
    all_frames_data_np = np.transpose(all_frames_data_np, (1, 2, 0))
    with open(filename, mode='w', newline='') as file:
        csv_writer = csv.writer(file, delimiter=',')
        for frame_index in range(all_frames_data_np.shape[2]):
            for point_index in range(all_frames_data_np.shape[0]):
                x, y, z = all_frames_data_np[point_index, :, frame_index]
                csv_writer.writerow([x,y,z])
            csv_writer.writerow([])

def extract_pose_from_video(video_path, csv_path):
    cap = cv2.VideoCapture(video_path)
    with mp_holistic.Holistic(smooth_landmarks=True, min_detection_confidence=0.5, min_tracking_confidence=0.5) as holistic:
        data = []
        while cap.isOpened():
            success, image = cap.read()
            if not success:
                break
            image.flags.writeable = False
            image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
            frame_data = extract_pose_landmarks(image, holistic)
            if len(frame_data) == 33:
                data.append(frame_data)
    cap.release()
    cv2.destroyAllWindows()
    save_landmarks_to_csv(data, csv_path)

  
#if __name__=='__main__':
    #extract_pose_from_video('static\uploads\front.avi','ViconMarkerless2Opensim\uploads\result.csv')
        