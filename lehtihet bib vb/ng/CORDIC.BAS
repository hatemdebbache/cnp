Attribute VB_Name = "CORDIC"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << CORDIC >> version 01
'  Requis   :
'  Externes :
'
'  Descrip. : Coordinate Rotational Digital Computation
'
' ================================================
'
'                  sinus et cosinus
'
' ================================================
Option Explicit
Private Const mPI2 As Double = 1.570796326795
Private mNPTS As Integer
Private mTAN() As Double
Private mPRODCOS As Double
Sub Cordic_Init(ByVal N As Integer)

    Dim phi As Double
    Dim k As Integer
    
    If N < 1 Then Stop
    mNPTS = N
    ReDim mTAN(N)
    
    mPRODCOS = 1#
    phi = mPI2
    For k = 1 To mNPTS
        phi = 0.5 * phi
        mTAN(k) = Tan(phi)
        mPRODCOS = mPRODCOS * Cos(phi)
    Next k
End Sub
Sub SinCos(ByVal angle As Double, _
           ByRef v As Double, _
           ByRef u As Double, _
  Optional ByRef r As Double)
    
    Dim k As Integer
    Dim phi As Double, ang As Double
    Dim c As Double, t As Double
    Dim temp As Double, dum As Double
    
    If Abs(angle) > mPI2 Then Stop
    
    ' Recherche des poids wi et cumul de A*
    phi = mPI2
    r = angle
    ang = 0#
    u = mPRODCOS
    v = 0#
    For k = 1 To mNPTS
        phi = 0.5 * phi
        If r > 0 Then
            r = r - phi
            ang = ang + phi
            dum = mTAN(k)
        Else
            r = r + phi
            ang = ang - phi
            dum = -mTAN(k)
        End If
        temp = u
        u = u - v * dum
        v = v + temp * dum
    Next k
    r = angle - ang
    dum = r * r
    c = 1# - 0.5 * dum
    dum = r * dum
    t = r * (1# + dum / 3#)
    temp = u
    u = c * (u - v * t)          ' cosA
    v = c * (v + temp * t)       ' sinA
End Sub

