package internal

import (
	"github.com/cloudflare/circl/sign/dilithium/internal/common"
)

// Writes p with norm less than or equal η into buf, which must be of
// size PolyLeqEtaSize.
//
// Assumes coefficients of p are not normalized, but in [q-η,q+η].
func PolyPackLeqEta(p *common.Poly, buf []byte) {
	if DoubleEtaBits == 4 { // compiler eliminates branch
		j := 0
		for i := 0; i < PolyLeqEtaSize; i++ {
			buf[i] = (byte(common.Q+Eta-p[j]) |
				byte(common.Q+Eta-p[j+1])<<4)
			j += 2
		}
	} else if DoubleEtaBits == 3 {
		j := 0
		for i := 0; i < PolyLeqEtaSize; i += 3 {
			buf[i] = (byte(common.Q+Eta-p[j]) |
				(byte(common.Q+Eta-p[j+1]) << 3) |
				(byte(common.Q+Eta-p[j+2]) << 6))
			buf[i+1] = ((byte(common.Q+Eta-p[j+2]) >> 2) |
				(byte(common.Q+Eta-p[j+3]) << 1) |
				(byte(common.Q+Eta-p[j+4]) << 4) |
				(byte(common.Q+Eta-p[j+5]) << 7))
			buf[i+2] = ((byte(common.Q+Eta-p[j+5]) >> 1) |
				(byte(common.Q+Eta-p[j+6]) << 2) |
				(byte(common.Q+Eta-p[j+7]) << 5))
			j += 8
		}
	} else {
		panic("eta not supported")
	}
}

// Sets p to the polynomial of norm less than or equal η encoded in the
// given buffer of size PolyLeqEtaSize.
//
// Output coefficients of p are not normalized, but in [q-η,q+η] provided
// buf was created using PackLeqEta.
//
// Beware, for arbitrary buf the coefficients of p might en up in
// the interval [q-2^b,q+2^b] where b is the least b with η≤2^b.
func PolyUnpackLeqEta(p *common.Poly, buf []byte) {
	if DoubleEtaBits == 4 { // compiler eliminates branch
		j := 0
		for i := 0; i < PolyLeqEtaSize; i++ {
			p[j] = common.Q + Eta - uint32(buf[i]&15)
			p[j+1] = common.Q + Eta - uint32(buf[i]>>4)
			j += 2
		}
	} else if DoubleEtaBits == 3 {
		j := 0
		for i := 0; i < PolyLeqEtaSize; i += 3 {
			p[j] = common.Q + Eta - uint32(buf[i]&7)
			p[j+1] = common.Q + Eta - uint32((buf[i]>>3)&7)
			p[j+2] = common.Q + Eta - uint32((buf[i]>>6)|((buf[i+1]<<2)&7))
			p[j+3] = common.Q + Eta - uint32((buf[i+1]>>1)&7)
			p[j+4] = common.Q + Eta - uint32((buf[i+1]>>4)&7)
			p[j+5] = common.Q + Eta - uint32((buf[i+1]>>7)|((buf[i+2]<<1)&7))
			p[j+6] = common.Q + Eta - uint32((buf[i+2]>>2)&7)
			p[j+7] = common.Q + Eta - uint32((buf[i+2]>>5)&7)
			j += 8
		}
	} else {
		panic("eta not supported")
	}
}

// Writes v with coefficients in {0, 1} of which at most ω non-zero
// to buf, which must have length ω+k.
func (v *VecK) PackHint(buf []byte) {
	// The packed hint starts with the indices of the non-zero coefficients
	// For instance:
	//
	//    (x⁵⁶ + x¹⁰⁰, x²⁵⁵, 0, x² + x²³, x¹)
	//
	// Yields
	//
	//  56, 100, 255, 2, 23, 1
	//
	// Then we pad with zeroes until we have a list of ω items:
	// //  56, 100, 255, 2, 23, 1, 0, 0, ..., 0
	//
	// Then we finish with a list of the switch-over-indices in this
	// list between polynomials, so:
	//
	//  56, 100, 255, 2, 23, 1, 0, 0, ..., 0, 2, 3, 3, 5, 6

	off := uint8(0)
	for i := 0; i < K; i++ {
		for j := uint16(0); j < common.N; j++ {
			if v[i][j] != 0 {
				buf[off] = uint8(j)
				off++
			}
		}
		buf[Omega+i] = off
	}
	for ; off < Omega; off++ {
		buf[off] = 0
	}
}

// Sets v to the vector encoded using VecK.PackHint()
//
// Returns whether unpacking was successful.
func (v *VecK) UnpackHint(buf []byte) bool {
	// A priori, there would be several reasonable ways to encode the same
	// hint vector.  We take care to only allow only one encoding, to ensure
	// "strong unforgeability".
	//
	// See PackHint() source for description of the encoding.
	*v = VecK{}         // zero v
	prevSOP := uint8(0) // previous switch-over-point
	for i := 0; i < K; i++ {
		SOP := buf[Omega+i]
		if SOP < prevSOP || SOP > Omega {
			return false // ensures switch-over-points are increasing
		}
		for j := prevSOP; j < SOP; j++ {
			if j > prevSOP && buf[j] <= buf[j-1] {
				return false // ensures indices are increasing (within a poly)
			}
			v[i][buf[j]] = 1
		}
		prevSOP = SOP
	}
	for j := prevSOP; j < Omega; j++ {
		if buf[j] != 0 {
			return false // ensures padding indices are zero
		}
	}

	return true
}

// Writes p whose coefficients are in norm less than γ₁ into buf
// which has to be of length PolyLeGamma1Size.
//
// Assumes p is normalized.
func PolyPackLeGamma1(p *common.Poly, buf []byte) {
	j := 0
	for i := 0; i < PolyLeGamma1Size; i += 5 {
		// Coefficients are in [0, γ₁) ∪ (Q-γ₁, Q)
		p0 := Gamma1 - 1 - p[j]         // ... in [0, γ₁) ∪ [γ₁-1-Q, 2(γ₁-1)-Q]
		p0 += uint32(int32(p0)>>31) & Q // ... in [0, 2(γ₁-1)]
		p1 := Gamma1 - 1 - p[j+1]
		p1 += uint32(int32(p1)>>31) & Q

		buf[i] = byte(p0)
		buf[i+1] = byte(p0 >> 8)
		buf[i+2] = byte(p0>>16) | byte(p1<<4)
		buf[i+3] = byte(p1 >> 4)
		buf[i+4] = byte(p1 >> 12)
		j += 2
	}
}

// Sets p to the polynomial packed into buf by PackLeGamma1.
//
// p will be normalized.
//
// Beware, for arbitrary buf the coefficients of p might exceed γ₁.
func PolyUnpackLeGamma1(p *common.Poly, buf []byte) {
	j := 0
	for i := 0; i < PolyLeGamma1Size; i += 40 {
		a0 := binary.LittleEndian.Uint64(buf[i:])
		a1 := binary.LittleEndian.Uint64(buf[i+8:])
		a2 := binary.LittleEndian.Uint64(buf[i+16:])
		a3 := binary.LittleEndian.Uint64(buf[i+24:])
		a4 := binary.LittleEndian.Uint64(buf[i+32:])

		p0 := Gamma1 - 1 - uint32(a0&0xfffff)
		p1 := Gamma1 - 1 - uint32((a0>>20)&0xfffff)
		p2 := Gamma1 - 1 - uint32((a0>>40)&0xfffff)
		p3 := Gamma1 - 1 - uint32(((a0>>60)|(a1<<4))&0xfffff)
		p4 := Gamma1 - 1 - uint32((a1>>16)&0xfffff)
		p5 := Gamma1 - 1 - uint32((a1>>36)&0xfffff)
		p6 := Gamma1 - 1 - uint32(((a1>>56)|(a2<<8))&0xfffff)
		p7 := Gamma1 - 1 - uint32((a2>>12)&0xfffff)
		p8 := Gamma1 - 1 - uint32((a2>>32)&0xfffff)
		p9 := Gamma1 - 1 - uint32(((a2>>52)|(a3<<12))&0xfffff)
		p10 := Gamma1 - 1 - uint32((a3>>8)&0xfffff)
		p11 := Gamma1 - 1 - uint32((a3>>28)&0xfffff)
		p12 := Gamma1 - 1 - uint32(((a3>>48)|(a4<<16))&0xfffff)
		p13 := Gamma1 - 1 - uint32((a4>>4)&0xfffff)
		p14 := Gamma1 - 1 - uint32((a4>>24)&0xfffff)
		p15 := Gamma1 - 1 - uint32((a4>>44)&0xfffff)

		p0 += uint32(int32(p0)>>31) & Q
		p1 += uint32(int32(p1)>>31) & Q
		p2 += uint32(int32(p2)>>31) & Q
		p3 += uint32(int32(p3)>>31) & Q
		p4 += uint32(int32(p4)>>31) & Q
		p5 += uint32(int32(p5)>>31) & Q
		p6 += uint32(int32(p6)>>31) & Q
		p7 += uint32(int32(p7)>>31) & Q
		p8 += uint32(int32(p8)>>31) & Q
		p9 += uint32(int32(p9)>>31) & Q
		p10 += uint32(int32(p10)>>31) & Q
		p11 += uint32(int32(p11)>>31) & Q
		p12 += uint32(int32(p12)>>31) & Q
		p13 += uint32(int32(p13)>>31) & Q
		p14 += uint32(int32(p14)>>31) & Q
		p15 += uint32(int32(p15)>>31) & Q

		p[j] = p0
		p[j+1] = p1
		p[j+2] = p2
		p[j+3] = p3
		p[j+4] = p4
		p[j+5] = p5
		p[j+6] = p6
		p[j+7] = p7
		p[j+8] = p8
		p[j+9] = p9
		p[j+10] = p10
		p[j+11] = p11
		p[j+12] = p12
		p[j+13] = p13
		p[j+14] = p14
		p[j+15] = p15

		j += 16
	}
}

