"""Replay Vim's alternate screen and ignore its nonvisual terminal probes."""
import copy
import pyte


class TerminalScreen(pyte.Screen):
    saved_main = None

    def set_mode(self, *modes, **kwargs):
        if kwargs.get("private") and 1049 in modes and self.saved_main is None:
            self.saved_main = copy.deepcopy((self.buffer, self.cursor, self.margins))
            self.erase_in_display(2)
            self.cursor_position()
        super().set_mode(*modes, **kwargs)

    def reset_mode(self, *modes, **kwargs):
        super().reset_mode(*modes, **kwargs)
        if kwargs.get("private") and 1049 in modes and self.saved_main is not None:
            self.buffer, self.cursor, self.margins = self.saved_main
            self.saved_main = None
            self.dirty.update(range(self.lines))

    def select_graphic_rendition(self, *attrs, private=False):
        if not private:
            super().select_graphic_rendition(*attrs)
