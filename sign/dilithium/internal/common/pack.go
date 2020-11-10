package common

// Sets p to the polynomial whose coefficients are less than 512 encoded
// into buf (which must be of size PolyT1Size).
//
// p will be normalized.
func (p *Poly) UnpackT1(buf []byte) {
	j := 0
	for i := 0; i < PolyT1Size; i += 9 {
		p[j] = (uint32(buf[i]) | (uint32(buf[i+1]) << 8)) & 0x1ff
		p[j+1] = (uint32(buf[i+1]>>1) | (uint32(buf[i+2]) << 7)) & 0x1ff
		p[j+2] = (uint32(buf[i+2]>>2) | (uint32(buf[i+3]) << 6)) & 0x1ff
		p[j+3] = (uint32(buf[i+3]>>3) | (uint32(buf[i+4]) << 5)) & 0x1ff
		p[j+4] = (uint32(buf[i+4]>>4) | (uint32(buf[i+5]) << 4)) & 0x1ff
		p[j+5] = (uint32(buf[i+5]>>5) | (uint32(buf[i+6]) << 3)) & 0x1ff
		p[j+6] = (uint32(buf[i+6]>>6) | (uint32(buf[i+7]) << 2)) & 0x1ff
		p[j+7] = (uint32(buf[i+7]>>7) | (uint32(buf[i+8]) << 1)) & 0x1ff
		j += 8
	}
}

// Writes p whose coefficients are in (-2ᵈ⁻¹, 2ᵈ⁻¹] into buf which
// has to be of length at least PolyT0Size.
//
// Assumes that the coefficients are not normalized, but lie in the
// range (q-2ᵈ⁻¹, q+2ᵈ⁻¹].
func (p *Poly) PackT0(buf []byte) {
	j := 0
	for i := 0; i < PolyT0Size; i += 7 {
		p0 := Q + (1 << (D - 1)) - p[j]
		p1 := Q + (1 << (D - 1)) - p[j+1]
		p2 := Q + (1 << (D - 1)) - p[j+2]
		p3 := Q + (1 << (D - 1)) - p[j+3]

		buf[i] = byte(p0)
		buf[i+1] = byte(p0>>8) | byte(p1<<6)
		buf[i+2] = byte(p1 >> 2)
		buf[i+3] = byte(p1>>10) | byte(p2<<4)
		buf[i+4] = byte(p2 >> 4)
		buf[i+5] = byte(p2>>12) | byte(p3<<2)
		buf[i+6] = byte(p3 >> 6)
		j += 4
	}
}

// Sets p to the polynomial packed into buf by PackT0.
//
// The coefficients of p will not be normalized, but will lie
// in (-2ᵈ⁻¹, 2ᵈ⁻¹].
func (p *Poly) UnpackT0(buf []byte) {
	j := 0
	for i := 0; i < PolyT0Size; i += 7 {
		p[j] = Q + (1 << (D - 1)) - ((uint32(buf[i]) |
			(uint32(buf[i+1]) << 8)) & 0x3fff)
		p[j+1] = Q + (1 << (D - 1)) - (((uint32(buf[i+1]) >> 6) |
			(uint32(buf[i+2]) << 2) |
			(uint32(buf[i+3]) << 10)) & 0x3fff)
		p[j+2] = Q + (1 << (D - 1)) - (((uint32(buf[i+3]) >> 4) |
			(uint32(buf[i+4]) << 4) |
			(uint32(buf[i+5]) << 12)) & 0x3fff)
		p[j+3] = Q + (1 << (D - 1)) - ((uint32(buf[i+5]) >> 2) |
			(uint32(buf[i+6]) << 6))
		j += 4
	}
}

// Writes p whose coefficients are less than 512 into buf, which must be
// of size at least PolyT1Size .
//
// Assumes coefficients of p are normalized.
func (p *Poly) PackT1(buf []byte) {
	j := 0
	for i := 0; i < PolyT1Size; i += 9 {
		buf[i] = byte(p[j])
		buf[i+1] = byte(p[j]>>8) | byte(p[j+1]<<1)
		buf[i+2] = byte(p[j+1]>>7) | byte(p[j+2]<<2)
		buf[i+3] = byte(p[j+2]>>6) | byte(p[j+3]<<3)
		buf[i+4] = byte(p[j+3]>>5) | byte(p[j+4]<<4)
		buf[i+5] = byte(p[j+4]>>4) | byte(p[j+5]<<5)
		buf[i+6] = byte(p[j+5]>>3) | byte(p[j+6]<<6)
		buf[i+7] = byte(p[j+6]>>2) | byte(p[j+7]<<7)
		buf[i+8] = byte(p[j+7] >> 1)
		j += 8
	}
}

// Writes p whose coefficients are in [0, 16) to buf, which must be of
// length N/2.
func (p *Poly) packLe16Generic(buf []byte) {
	j := 0
	for i := 0; i < PolyLe16Size; i++ {
		buf[i] = byte(p[j]) | byte(p[j+1]<<4)
		j += 2
	}
}

// Writes p with 60 non-zero coefficients {-1,1} to buf, which must have
// length 40.
func (p *Poly) PackB60(buf []byte) {
	// We start with a mask of the non-zero positions of p (which is 32 bytes)
	// and then append 60 packed bits, where a one indicates a negative
	// coefficients.
	var signs uint64
	mask := uint64(1)
	for i := 0; i < 32; i++ {
		buf[i] = 0
		for j := 0; j < 8; j++ {
			if p[8*i+j] != 0 {
				buf[i] |= 1 << uint(j)
				if p[8*i+j] == Q-1 {
					signs |= mask
				}
				mask <<= 1
			}
		}
	}
	for i := uint64(0); i < 8; i++ {
		buf[i+32] = uint8(signs >> (8 * i))
	}
}

// UnpackB60 sets p to the polynomial packed into buf with Poly.PackB60().
//
// Returns whether unpacking was successful.
func (p *Poly) UnpackB60(buf []byte) bool {
	*p = Poly{} // zero p
	signs := (uint64(buf[32]) | (uint64(buf[33]) << 8) |
		(uint64(buf[34]) << 16) | (uint64(buf[35]) << 24) |
		(uint64(buf[36]) << 32) | (uint64(buf[37]) << 40) |
		(uint64(buf[38]) << 48) | (uint64(buf[39]) << 56))
	if signs>>60 != 0 {
		return false // ensure unused bits are zero for strong unforgeability
	}

	for i := 0; i < 32; i++ {
		for j := 0; j < 8; j++ {
			if (buf[i]>>uint(j))&1 == 1 {
				p[8*i+j] = 1
				// Note 1 ^ (1 | (Q-1)) = Q-1 and (-1)&x = x
				p[8*i+j] ^= uint32(-(signs & 1)) & (1 | (Q - 1))
				signs >>= 1
			}
		}
	}

	return true
}
