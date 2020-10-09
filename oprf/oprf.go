package oprf

import (
	"encoding/binary"
	"errors"
	"fmt"

	"github.com/cloudflare/circl/oprf/group"
)

const (
	// OPRFCurve25519 is the constant to represent the OPRF curve25519 with SHA-512 (ELL2-RO) group.
	OPRFCurve25519 uint16 = 0x0001
	// OPRFCurve448 is the constant to represent the OPRF curve448 with SHA-512 (ELL2-RO) group.
	OPRFCurve448 uint16 = 0x0002
	// OPRFP256 is the constant to represent the OPRF P-256 with SHA-512 (SSWU-RO) group.
	OPRFP256 uint16 = 0x0003
	// OPRFP384 is the constant to represent the OPRF P-384 with SHA-512 (SSWU-RO) group.
	OPRFP384 uint16 = 0x0004
	// OPRFP521 is the constant to represent the OPRF P-521 with SHA-512 (SSWU-RO) group.
	OPRFP521 uint16 = 0x0005
)

const (
	// OPRFMode is the context string to define a OPRF.
	OPRFMode string = "000"
)

var (
	// ErrUnsupportedGroup is an error stating that the ciphersuite chosen is not supported
	ErrUnsupportedGroup = errors.New("the chosen group is not supported")
)

// BlindToken corresponds to a token that has been blinded.
// Internally, it is a serialized Element.
type BlindToken []byte

// IssuedToken corresponds to a token that has been issued.
// Internally, it is a serialized Element.
type IssuedToken []byte

// Token is the object issuance of the protocol.
type Token struct {
	data  []byte
	blind *group.Scalar
}

// Evaluation corresponds to the evaluation over a token.
type Evaluation struct {
	element []byte
}

// KeyPair is an struct containing a public and private key.
type KeyPair struct {
	PubK  *group.Element
	PrivK *group.Scalar
}

// Client is a representation of a Client during protocol execution.
type Client struct {
	suite *group.Ciphersuite
	ctx   string
}

// ClientContext implements the functionality of a Client.
type ClientContext interface {
	// Blind generates a token and blinded data.
	Blind(in []byte) (*Token, *BlindToken, error)

	// Unblind unblinds the server response.
	Unblind(t *Token, e *Evaluation) (IssuedToken, error)

	// Finalize outputs a byte array that corresponds to its input.
	Finalize(t *Token, issuedT IssuedToken, info []byte) []byte
}

// Server is a representation of a Server during protocol execution.
type Server struct {
	suite *group.Ciphersuite
	ctx   string
	K     *KeyPair
}

// ServerContext implements the functionality of a Server.
// TODO: add FullEvaluate
type ServerContext interface {
	// Evaluate evaluates the token.
	Evaluate(b BlindToken) *Evaluation
}

func generateCtx(suiteID uint16) string {
	ctx := OPRFMode + fmt.Sprintf("%x", suiteID)

	return ctx
}

func generateKeys(suite *group.Ciphersuite) (*KeyPair, error) {
	privK, err := suite.RandomScalar()
	if err != nil {
		return nil, err
	}

	pubK := suite.ScalarMultBase(privK)

	return &KeyPair{pubK, privK}, nil
}

func suiteFromID(suiteID uint16) (*group.Ciphersuite, error) {
	var err error
	var suite *group.Ciphersuite

	// TODO: add other suites
	switch suiteID {
	case OPRFP256:
		suite, err = group.NewSuite("P-256")
	case OPRFP384:
		suite, err = group.NewSuite("P-384")
	case OPRFP521:
		suite, err = group.NewSuite("P-521")
	default:
		return suite, ErrUnsupportedGroup
	}
	if err != nil {
		return nil, err
	}

	return suite, nil
}

// NewServer creates a new instantiation of a Server.
func NewServer(suiteID uint16) (*Server, error) {
	suite, err := suiteFromID(suiteID)
	if err != nil {
		return nil, err
	}

	ctx := generateCtx(suiteID)

	keyPair, err := generateKeys(suite)
	if err != nil {
		return nil, err
	}

	server := &Server{}
	server.suite = suite
	server.ctx = ctx
	server.K = keyPair

	return server, nil
}

// Evaluate creates an evaluation of the blided token.
func (s *Server) Evaluate(b BlindToken) *Evaluation {
	p := group.NewElement(s.suite)
	p.Deserialize(b)

	z := p.ScalarMult(s.K.PrivK)

	ser := z.Serialize()

	eval := &Evaluation{ser}
	return eval
}

// NewClient creates a new instantiation of a Client.
func NewClient(suiteID uint16) (*Client, error) {
	suite, err := suiteFromID(suiteID)
	if err != nil {
		return nil, err
	}

	ctx := generateCtx(suiteID)

	client := &Client{}
	client.suite = suite
	client.ctx = ctx

	return client, nil
}

// Blind generates a token and blinded data.
func (c *Client) Blind(in []byte) (*Token, BlindToken, error) {
	r, err := c.suite.RandomScalar()
	if err != nil {
		return nil, nil, errors.New("failed at blinding")
	}

	p, err := c.suite.HashToGroup(in)
	if err != nil {
		return nil, nil, errors.New("failed at blinding")
	}

	t := p.ScalarMult(r)
	bToken := t.Serialize()

	token := &Token{in, r}
	return token, bToken, nil
}

// Unblind unblinds the server response.
func (c *Client) Unblind(t *Token, e *Evaluation) (IssuedToken, error) {
	p := group.NewElement(c.suite)
	p.Deserialize(e.element)

	r := t.blind
	rInv := r.Inv()

	tt := p.ScalarMult(rInv)
	iToken := tt.Serialize()

	return IssuedToken(iToken), nil
}

// Finalize outputs a byte array that corresponds to the client input.
func (c *Client) Finalize(t *Token, issuedT IssuedToken, info []byte) []byte {
	h := c.suite.Hash

	lenBuf := make([]byte, 2)
	binary.BigEndian.PutUint16(lenBuf, uint16(len(t.data)))
	_, _ = h.Write(lenBuf)
	_, _ = h.Write(t.data)

	binary.BigEndian.PutUint16(lenBuf, uint16(len(issuedT)))
	_, _ = h.Write(lenBuf)
	_, _ = h.Write(issuedT)

	binary.BigEndian.PutUint16(lenBuf, uint16(len(info)))
	_, _ = h.Write(lenBuf)
	_, _ = h.Write(info)

	dst := []byte("RFCXXXX-Finalize")
	dst = append(dst, c.ctx...)
	binary.BigEndian.PutUint16(lenBuf, uint16(len(dst)))
	_, _ = h.Write(lenBuf)
	_, _ = h.Write(dst)

	return h.Sum(nil)
}
