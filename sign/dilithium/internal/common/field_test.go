package common

import (
	"math/rand"
	"testing"
)

func TestModQ(t *testing.T) {
	for i := 0; i < 1000; i++ {
		x := rand.Uint32()
		y := modQ(x)
		if y > Q {
			t.Fatalf("modQ(%d) > Q", x)
		}
		if y != x%Q {
			t.Fatalf("modQ(%d) != %d (mod Q)", x, x)
		}
	}
}

func TestReduceLe2Q(t *testing.T) {
	for i := 0; i < 1000; i++ {
		x := rand.Uint32()
		y := reduceLe2Q(x)
		if y > 2*Q {
			t.Fatalf("reduce_le2q(%d) > 2Q", x)
		}
		if y%Q != x%Q {
			t.Fatalf("reduce_le2q(%d) != %d (mod Q)", x, x)
		}
	}
}

func TestPower2Round(t *testing.T) {
	for a := uint32(0); a < Q; a++ {
		a0PlusQ, a1 := power2round(a)
		a0 := int32(a0PlusQ) - int32(Q)
		if int32(a) != a0+int32((1<<D)*a1) {
			t.Fatalf("power2round(%v) doesn't recombine", a)
		}
		if (-(1 << (D - 1)) >= a0) || (a0 > 1<<(D-1)) {
			t.Fatalf("power2round(%v): a0 out of bounds", a)
		}
		if a1 > (1 << (23 - 14)) { // 23 ~ 2log Q.
			t.Fatalf("power2round(%v): a1 out of bounds", a)
		}
	}
}
