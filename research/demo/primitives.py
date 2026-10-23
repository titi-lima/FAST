"""A compact, data-driven shingles/MinHash illustration using text and rectangles."""
import math
from PIL import Image, ImageDraw, ImageFont

BG = "#171717"
FG = "#dedede"
MUTED = "#92949d"
BLUE = "#8caaee"
GREEN = "#a6d189"
ORANGE = "#e5a67a"
COLORS = {"A": BLUE, "B": ORANGE, "C": GREEN}
DURATION = 26


class Illustration:
    def __init__(self, data, font):
        self.data = data
        self.font = ImageFont.truetype(font, 26)
        self.small = ImageFont.truetype(font, 21)
        self.tiny = ImageFont.truetype(font, 18)
        self.title = ImageFont.truetype(font, 30)
        self.cw = self.font.getlength("M")

    def frame(self, seconds):
        im = Image.new("RGB", (1920, 1080), BG)
        d = ImageDraw.Draw(im)
        data = self.data
        stage = 0 if seconds < 6 else 1 if seconds < 14 else 2
        heading = ["Test names → 5-character shingles", "MinHash → one minimum per hash function",
                   "FAST-pw · compare candidates"][stage]
        d.text((80, 54), heading, font=self.title, fill=FG)
        if stage == 2:
            subtitle = "A selected · lower estimated similarity next" if seconds < 19 else "A + C selected · compare the remaining test"
            d.text((80, 108), subtitle, font=self.small, fill=MUTED)
        for row, key in enumerate(("A", "B", "C")):
            y = 190 + row * 220
            test = data["tests"][key]
            text = test["nodeid"]
            d.text((80, y), key, font=self.font, fill=COLORS[key])
            d.text((140, y), text, font=self.font, fill=FG)
            if stage == 0:
                start = text.index("shipping") if key != "C" else text.index("test_tax")
                offset = min(3, int(max(0, seconds - 1) / 1.1))
                index = start + offset
                x = 140 + index * self.cw
                d.rectangle((x - 3, y - 3, x + 5 * self.cw + 2, y + 35), outline=COLORS[key], width=2)
                for chip in range(min(4, offset + 1)):
                    shingle = text[start + chip:start + chip + 5]
                    cx = 140 + chip * 180
                    d.rectangle((cx, y + 65, cx + 152, y + 115), outline=COLORS[key], width=1)
                    d.text((cx + 26, y + 72), shingle, font=self.font, fill=FG)
                d.text((875, y + 76), "…", font=self.font, fill=MUTED)
                continue
            # Actual signature values are abbreviated as four hexadecimal digits;
            # colors and equality counts always use the complete integer values.
            ready = 1 if stage == 2 else min(1, max(0, (seconds - 9 - row * .8) / 1.8))
            for col, value in enumerate(test["hashvalues"]):
                if col >= math.ceil(10 * ready):
                    continue
                x = 1060 + col * 72
                reference = data["cumulative_ac"] if seconds >= 22 and key == "B" else data["tests"]["A"]["hashvalues"]
                a = reference[col]
                stroke = BLUE if key == "A" or value == a else COLORS[key]
                compared = 10
                if stage == 2 and key != "A":
                    compared = min(10, max(0, int((seconds - (22 if seconds >= 22 and key == "B" else 14)) * (5 if seconds >= 22 and key == "B" else 3.5))))
                    if col >= compared:
                        stroke = MUTED
                fill = "#253145" if stage == 2 and key != "A" and value == a and col < compared else BG
                d.rectangle((x, y + 64, x + 64, y + 121), outline=stroke, fill=fill, width=2)
                d.text((x + 9, y + 79), f"{value:08x}"[:4], font=self.tiny, fill=stroke)
                if row == 0:
                    d.text((x + 13, y + 128), f"h{col + 1}", font=self.tiny, fill=MUTED)
            if stage == 1 and key == "A":
                for index, shingle in enumerate(data["shingle_examples"]):
                    x = 140 + index * 180
                    stroke = GREEN if index == 0 and seconds >= 8 else MUTED
                    d.rectangle((x, y + 65, x + 152, y + 115), outline=stroke, width=2 if stroke == GREEN else 1)
                    d.text((x + 26, y + 72), shingle["text"], font=self.font, fill=FG)
                    if seconds >= 7:
                        d.text((x + 19, y + 128), f"{shingle['h1']:08x}", font=self.small, fill=stroke)
                d.text((875, y + 76), "…", font=self.font, fill=MUTED)
                if seconds >= 8:
                    d.text((140, y + 172), "h1 · minimum", font=self.small, fill=GREEN)
                    d.line((220, y + 210, 1090, y + 210, 1090, y + 122), fill=GREEN, width=2)
                    d.line((1083, y + 134, 1090, y + 122, 1097, y + 134), fill=GREEN, width=2)
            if stage == 2 and key != "A":
                reference = data["cumulative_ac"] if seconds >= 22 and key == "B" else data["tests"]["A"]["hashvalues"]
                matches = sum(a == b for a, b in zip(test["hashvalues"][:compared], reference[:compared]))
                label = "selected" if key == "C" and seconds >= 18 else (f"{matches}/{compared} matching minima" if compared else "compare minima")
                d.text((140, y + 73), label, font=self.small, fill=COLORS[key])
        if stage == 1:
            d.text((1060, 876), "10 minima · hex prefixes", font=self.small, fill=MUTED)
        if seconds >= 19:
            d.text((140, 814), "merged minima · A + C", font=self.small, fill=MUTED)
            for col in range(min(10, int((seconds - 19) * 6))):
                value = data["cumulative_ac"][col]
                x = 1060 + col * 72
                stroke = BLUE if value == data["tests"]["A"]["hashvalues"][col] else GREEN
                d.rectangle((x, 800, x + 64, 857), outline=stroke, width=2)
                d.text((x + 9, 815), f"{value:08x}"[:4], font=self.tiny, fill=stroke)
        if seconds >= 14:
            d.text((80, 922), "order", font=self.small, fill=MUTED)
            visible = 1 + int(seconds >= 18) + int(seconds >= 24.3)
            for index, key in enumerate(data["order"][:visible]):
                x = 300 + index * 420
                d.rectangle((x, 907, x + 90, 980), outline=COLORS[key], width=2)
                d.text((x + 36, 922), key, font=self.title, fill=COLORS[key])
                if index:
                    d.line((x - 310, 945, x - 18, 945), fill=MUTED, width=2)
                    d.line((x - 28, 938, x - 18, 945, x - 28, 952), fill=MUTED, width=2)
        return im
