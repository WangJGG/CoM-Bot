# coding:utf-8

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.responses import RedirectResponse
from fastapi.middleware.wsgi import WSGIMiddleware
from flask import Flask,render_template,request,redirect,url_for
from werkzeug.utils import secure_filename
import subprocess
import os

flask_app = Flask(__name__)


@flask_app.route('/')
def exercisechoice():
    return render_template('exercisechoice.html')


@flask_app.route('/exercise', methods=['POST'])
def exercise():
    exercise_name = request.form.get('exercise')
    # 根据选中的运动跳转到对应的训练页面
    if exercise_name == '下蹲':
        return redirect(url_for('squat'))
    elif exercise_name == '举重':
        return redirect(url_for('weightlifting'))
    elif exercise_name == '瑜伽':
        return redirect(url_for('yoga'))
    else:
        return '未知运动'


@flask_app.route('/squat')
def squat():
    return render_template('squat.html')


@flask_app.route('/weightlifting')
def weightlifting():
    return render_template('weightlifting.html')


@flask_app.route('/yoga')
def yoga():
    return render_template('yoga.html')

@flask_app.route('/upload', methods=['POST', 'GET'])  #主页面
def upload():
    if request.method == 'POST':
        f = request.files['file']
        basepath = os.path.dirname(__file__)  
        upload_path = os.path.join(basepath, 'static/uploads',secure_filename(f.filename)) 
        f.save(upload_path)
        return redirect(url_for('upload'))
    return render_template('upload.html')

if __name__ == '__main__':
    flask_app.run(debug=True)

app = FastAPI()

media_directory = "Result"
app.mount("/static", StaticFiles(directory=media_directory), name="static")

@app.get("/show-media", response_class=HTMLResponse)
async def show_media():
    try:
        media_files = os.listdir(media_directory)
        media_html = "<html><body><h3>反馈结果：</h3>"
        for media_file in media_files:
            if media_file.endswith((".mp4", ".avi", ".txt")):
                if media_file.endswith((".mp4", ".avi")):
                    media_html += f"<video src='/static/{media_file}' controls style='width:500px;'></video><br/>"
                else:
                    file_path = f'Result/{media_file}'
                    with open(file_path, "r") as txt_file:file_content = txt_file.read()
                    media_html += f"<p>{file_content}</p><br/>"
            else:
                media_html += f"<p>{media_file}</p><br/>"
        media_html += "</body></html>"
        return media_html
    except Exception as e:
        return HTMLResponse(content=f"<html><body><h3>Error Listing Media</h3><pre>{str(e)}</pre></body></html>", status_code=500)

@app.get("/run-matlab-and-PE")
async def run_matlab_and_PE():
    try:
        scripts = [
            (".\ViconMarkerless2Opensim\mediapipe_opensim.m", "Script 1"),
            (".\PiezoElectirc\PiezoElectric.m", "Script 2")]
        error_messages = []
        for script_path, script_name in scripts:
            command = f"matlab -batch \"run('{script_path}');exit;\""
            result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                error_messages.append(f"Error in {script_name}: {result.stderr}")
        if error_messages:
            return HTMLResponse(content=f"<html><body><h3>Failed to Execute MATLAB Script.</h3><pre>{'<br>'.join(error_messages)}</pre></body></html>", status_code=500)
        else:
         return RedirectResponse(url="/show-media", status_code=303)
    except Exception as e:
        return HTMLResponse(content=f"<html><body><h3>Error</h3><pre>{str(e)}</pre></body></html>", status_code=500)


@app.get("/", response_class=HTMLResponse)
def read_root():
    return """
    <html>
        <body>
            <h3>训练反馈生成中...</h3>
            <meta http-equiv="refresh" content="1; URL='/run-matlab-and-PE'">
        </body>
    </html>
    """

# 将Flask应用作为子应用挂载到FastAPI
app.mount("/upload", WSGIMiddleware(flask_app))