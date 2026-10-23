"""Render genuine terminal actions and computed primitives to exactly 3:00."""
import json
import math
import os
from pathlib import Path
import subprocess
import sys

from primitives import Illustration, DURATION as VISUAL_DURATION
from terminal import TerminalScreen
import pyte
from PIL import Image, ImageDraw, ImageFont

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(os.environ.get("FAST_DEMO_OUTPUT",ROOT/"output/video/three-minute"))
FPS=20
COLORS={"default":"#dedede","black":"#171717","red":"#ef767a",
        "green":"#a6d189","brown":"#e5c890","yellow":"#e5c890",
        "blue":"#8caaee","magenta":"#ca9ee6","cyan":"#81c8be","white":"#dedede",
        "brightblack":"#737994","brightred":"#ef9a9d","brightgreen":"#b5e3a4",
        "brightbrown":"#f2d7aa","brightyellow":"#f2d7aa","brightblue":"#a8bdf4",
        "brightmagenta":"#d7b6ed","brightcyan":"#a0ddd6","brightwhite":"#ffffff"}


def color(value,default):
    if value=="default":return default
    if value in COLORS:return COLORS[value]
    if len(value)==6:
        try:int(value,16);return "#"+value
        except ValueError:pass
    return default


def main():
    e=json.loads((OUT/"evidence.json").read_text())
    assert e["verified"]
    records=[json.loads(line) for line in (OUT/"fast-tcp-demo.cast").read_text().splitlines()]
    header,events=records[0],records[1:]
    target_duration = 180
    terminal_duration = target_duration - VISUAL_DURATION
    source_duration = e["duration_s"]
    pacing=json.loads((OUT/"pacing.json").read_text())["intervals"]
    segments=[]
    cursor=effective_duration=0
    def segment(start,end,speed):
        nonlocal effective_duration
        if end>start:
            length=(end-start)/speed
            segments.append((start,end,speed,effective_duration,effective_duration+length))
            effective_duration+=length
    for item in sorted(pacing,key=lambda item:item["start"]):
        assert cursor<=item["start"]<item["end"]<=source_duration and item["speed"]>=1
        segment(cursor,item["start"],1)
        segment(item["start"],item["end"],item["speed"])
        cursor=item["end"]
    segment(cursor,source_duration,1)
    def effective_time(source):
        for start,end,speed,a,b in segments:
            if source<=end:return a+(source-start)/speed
        return effective_duration
    def source_time(effective):
        for start,end,speed,a,b in segments:
            if effective<=b:return start+(effective-a)*speed
        return source_duration
    for check in range(101):
        source=source_duration*check/100
        assert abs(source_time(effective_time(source))-source)<1e-6
    if effective_duration < terminal_duration:
        raise SystemExit("Capture is too short; add useful recorded actions before rendering")
    source_insertion = next(scene["start"] for scene in e["scenes"] if scene["scene"] == "add regression tests")
    insertion = effective_time(source_insertion) / effective_duration * terminal_duration
    visual_data = json.loads((OUT/"visual-data.json").read_text())
    assert visual_data["verified"]
    font_path = os.environ.get("FAST_DEMO_FONT", "/System/Library/Fonts/Menlo.ttc")
    illustration = Illustration(visual_data, font_path)
    screen=TerminalScreen(header["width"],header["height"])
    stream=pyte.Stream(screen)
    regular=ImageFont.truetype(font_path,26)
    try:
        bold=ImageFont.truetype(font_path,26,index=1)
    except OSError:
        bold=regular
    cw=regular.getlength("M")
    line_height=32
    x0=round((1920-header["width"]*cw)/2)
    y0=28
    output=OUT/"fast-tcp-demo.mp4"
    encoder=subprocess.Popen(["ffmpeg","-hide_banner","-loglevel","warning","-y",
        "-f","rawvideo","-pix_fmt","rgb24","-s","1920x1080","-r",str(FPS),"-i","-",
        "-an","-c:v","libx264","-preset","veryfast","-crf","18","-pix_fmt","yuv420p",
        "-movflags","+faststart","-metadata","title=fast-tcp — three-minute demonstration",
        "-metadata","comment=Genuine shell/Vim capture, with a code-entry time lapse; a computed shingles/MinHash illustration uses text and rectangles. Constructed fixture; not a benchmark.",str(output)],stdin=subprocess.PIPE)
    index=0
    cached=None
    def movie_time(source):
        result = effective_time(source) / effective_duration * terminal_duration
        return result + VISUAL_DURATION if result >= insertion else result
    checkpoints={round(movie_time((e["scenes"][i+1]["start"] if i+1<len(e["scenes"]) else source_duration)-.7)*FPS):s["scene"].replace(" ","-") for i,s in enumerate(e["scenes"])}
    for frame in range(target_duration*FPS):
        seconds=frame/FPS
        if insertion <= seconds < insertion + VISUAL_DURATION:
            image=illustration.frame(seconds-insertion)
            encoder.stdin.write(image.tobytes())
            for moment in (4, 12, 19, 25):
                if abs(seconds-insertion-moment)<.5/FPS:
                    image.save(OUT/f"review-visual-{moment}.png")
            continue
        source_seconds=source_time((seconds if seconds<insertion else seconds-VISUAL_DURATION) / terminal_duration * effective_duration)
        changed=False
        while index<len(events) and events[index][0]<=source_seconds:
            _,kind,data=events[index]
            if kind=="o":stream.feed(data);changed=True
            index+=1
        if changed or cached is None:
            cached=Image.new("RGB",(1920,1080),"#171717")
            d=ImageDraw.Draw(cached)
            for row in range(header["height"]):
                col=0
                while col<header["width"]:
                    cell=screen.buffer[row][col]
                    style=(cell.fg,cell.bg,cell.bold,cell.reverse,cell.underscore)
                    end=col+1
                    while end<header["width"]:
                        next_cell=screen.buffer[row][end]
                        if (next_cell.fg,next_cell.bg,next_cell.bold,next_cell.reverse,next_cell.underscore)!=style:break
                        end+=1
                    chars="".join(screen.buffer[row][c].data for c in range(col,end))
                    fg=color(cell.fg,"#dedede")
                    bg=color(cell.bg,"#171717")
                    if cell.reverse:fg,bg=bg,fg
                    x,y=x0+col*cw,y0+row*line_height
                    if bg!="#171717":d.rectangle((x,y,x0+end*cw,y+line_height),fill=bg)
                    if chars.strip():d.text((x,y),chars,font=bold if cell.bold else regular,fill=fg)
                    if cell.underscore:d.line((x,y+29,x0+end*cw,y+29),fill=fg)
                    col=end
        image=cached.copy()
        if not screen.cursor.hidden and int(seconds*2)%2==0:
            d=ImageDraw.Draw(image)
            cx=x0+min(screen.cursor.x,header["width"]-1)*cw
            cy=y0+screen.cursor.y*line_height
            d.rectangle((cx,cy+2,cx+cw-1,cy+29),outline="#d4d4d4",width=2)
        encoder.stdin.write(image.tobytes())
        if frame in checkpoints:image.save(OUT/("review-"+checkpoints[frame]+".png"))
        if frame%(10*FPS)==0:print(f"RENDERED {seconds:.0f}s",flush=True)
    encoder.stdin.close()
    assert encoder.wait()==0
    manifest={"path":str(output),"duration_s":target_duration,"source_duration_s":source_duration,
              "terminal_playback_rate":effective_duration/terminal_duration,"code_entry_pacing":pacing,"visual_start_s":insertion,
              "visual_duration_s":VISUAL_DURATION,"fps":FPS,"resolution":[1920,1080],
              "cards":False,"clock_overlay":False,"fabricated_terminal_output":False}
    (OUT/"render-manifest.json").write_text(json.dumps(manifest,indent=2))
    print(json.dumps(manifest),flush=True)


if __name__=="__main__":main()
