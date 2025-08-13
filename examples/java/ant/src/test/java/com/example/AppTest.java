package com.example;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class AppTest {
    @Test
    void greetWithName() {
        String msg = App.greet("FAST");
        assertTrue(msg.contains("FAST"));
    }

    @Test
    void greetDefault() {
        String msg = App.greet("");
        assertEquals("Hello, World", msg);
    }
}


