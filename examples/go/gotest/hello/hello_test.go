package hello

import "testing"

func TestGreetWorld(t *testing.T) {
	got := Greet("")
	if got != "Hello, world!" {
		t.Fatalf("expected 'Hello, world!', got %q", got)
	}
}

func TestGreetName(t *testing.T) {
	got := Greet("Tiago")
	if got != "Hello, Tiago!" {
		t.Fatalf("expected 'Hello, Tiago!', got %q", got)
	}
}

func TestGreetAnother(t *testing.T) {
	got := Greet("FAST")
	if got != "Hello, FAST!" {
		t.Fatalf("expected 'Hello, FAST!', got %q", got)
	}
}
