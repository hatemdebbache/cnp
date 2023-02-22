Attribute VB_Name = "COEFF"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << COEFF >>
'  Requis   :          << SYSLIN >>
'  Externes :
'
'  Descrip. : Coefficients d'approximation
'
Option Explicit
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'      Coef. d'interpolation polynomiale simple
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Sub Newton(ByVal k1 As Integer, _
           ByVal k2 As Integer, _
           ByRef x() As Double, _
           ByRef y() As Double, _
           ByRef Tau() As Double)

    Dim k As Integer, m As Integer
    Dim ko As Integer, N As Integer
    Dim h As Double
    
    ko = k1 - 1
    N = k2 - ko
    If (N < 1) Then Stop ' Données non valides
    ReDim Tau(N)
    For k = 1 To N
        Tau(k) = y(k + ko)
    Next k
    If N = 1 Then Exit Sub
    For m = 1 To N - 1
        ko = k2 - N
        For k = N To m + 1 Step -1
            h = x(k + ko) - x(k + ko - m)
            If (h = 0) Then Stop ' Points non distincts
            Tau(k) = (Tau(k) - Tau(k - 1)) / h
        Next k
    Next m
End Sub
Sub Interpol(ByVal k1 As Integer, _
             ByVal k2 As Integer, _
             ByRef x() As Double, _
             ByRef y() As Double, _
             ByRef c() As Double)

    Dim r() As Double
    Dim k As Integer, j As Integer
    Dim ko As Integer, N As Integer, m As Integer
    Dim u As Double, v As Double, w As Double
    
    ko = k1 - 1          ' Indice de décalage
    m = k2 - k1          ' Degré du polynôme
    N = m + 1            ' Nombre total de points
    If (N < 1) Then Stop ' Données non valides
    ReDim c(N), r(m)
    w = 1#: u = x(k1)
    For k = 1 To m
        r(k) = -u * w
        w = 1#: u = x(k + k1)
        For j = 1 To k
            v = r(j)
            r(j) = v - u * w
            w = v
        Next j
    Next k
    For k = 1 To N
        c(k) = 0#
    Next k
    For k = 1 To N
        u = x(k + ko)
        w = N
        For j = 1 To m
            w = u * w + (N - j) * r(j)
        Next j
        v = y(k + ko) / w
        c(N) = c(N) + v
        w = 1#
        For j = 1 To m
            w = u * w + r(j)
            c(N - j) = c(N - j) + v * w
        Next j
    Next k
End Sub
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'  Coef. d'interpolation polynomiale par morceaux
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Sub Morceau(ByVal N As Integer, _
            ByRef x() As Double, _
            ByVal xo As Double, _
            ByRef k1 As Integer, _
            ByRef k2 As Integer, _
            ByRef h As Double, _
            ByRef v As Double, _
            ByRef w As Double, _
   Optional ByRef u1 As Double, _
   Optional ByRef u2 As Double)

    Dim m As Integer
    
    If (N < 2) Then Stop
    If (xo < x(2)) Then
        k1 = 1: k2 = 2
    ElseIf (xo >= x(N - 1)) Then
        k1 = N - 1: k2 = N
    Else
        k1 = 1: k2 = N
        Do While (k2 - k1) > 1
            m = (k1 + k2) / 2
            If (x(m) <= xo) Then k1 = m Else k2 = m
        Loop
    End If
    h = x(k2) - x(k1)
    If (h <= 0) Then Stop
    u1 = xo - x(k1): u2 = x(k2) - xo
    v = u1 / h:   w = 1# - v
End Sub
Sub Spline(ByVal MDF1 As Integer, _
           ByVal MDFn As Integer, _
           ByVal N As Integer, _
           ByRef x() As Double, _
           ByRef y() As Double, _
           ByRef s() As Double, _
  Optional ByVal d1 As Double = 0#, _
  Optional ByVal dn As Double = 0#)

    Dim z()  As Double
    Dim k    As Integer, bool As Boolean
    Dim h1   As Double, h2  As Double
    Dim dd1  As Double, ddn As Double
    Dim c As Double, b As Double, w As Double
                 
    If (N <= 2) Then Stop
    h1 = x(2) - x(1)
    If (h1 <= 0#) Then Stop
    ReDim z(1 To N)
    bool = ((MDF1 < 0) Or (MDFn < 0))
    
    If (MDF1 = 0) Then
        z(1) = 0#: s(1) = 0#
    Else
        If (MDF1 < 0) Then
            If (N < 5) Then
                dd1 = (-3# * y(1) + 4# * y(2) - y(3)) / (2# * h1)
            Else
                dd1 = (-25# * y(1) + 48# * y(2) - 36# * y(3) + 16# * y(4) - 3# * y(5)) / (12# * h1)
            End If
        Else
            dd1 = d1
        End If
        z(1) = 3# * ((y(2) - y(1)) / h1 - dd1) / h1
        s(1) = -0.5
    End If
    
    For k = 2 To N - 1
        h2 = x(k + 1) - x(k)
        If (h2 <= 0#) Then Stop
   '     If ((h1 <> h2) And (bool)) Then Stop
        c = h1 / (h1 + h2)
        w = 1# / (c * s(k - 1) + 2#)
        s(k) = w * (c - 1#)
        b = (c * y(k + 1) - y(k) + (1# - c) * y(k - 1)) * 6# / (h1 * h2)
        z(k) = (b - c * z(k - 1)) * w
        h1 = h2
    Next k

    If (MDFn = 0) Then
        s(N) = 0#: z(N) = 0#
    Else
        If (MDFn < 0) Then
            If (N < 5) Then
                ddn = (3# * y(N) - 4# * y(N - 1) + y(N - 2)) / (2# * h1)
            Else
                ddn = (25# * y(N) - 48# * y(N - 1) + 36# * y(N - 2) - 16# * y(N - 3) + 3# * y(N - 4)) / (12# * h1)
            End If
        Else
            ddn = dn
        End If
        z(N) = 6# * (ddn - (y(N) - y(N - 1)) / h1) / h1
        s(N) = (z(N) - z(N - 1)) / (s(N - 1) + 2#)
    End If

    For k = N - 1 To 1 Step -1
        s(k) = z(k) + s(k) * s(k + 1)
    Next k
    err = 0
End Sub
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'        Coef. d'approximation de moindres carrés
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Sub MC_Droite(ByVal k1 As Integer, _
              ByVal k2 As Integer, _
              ByRef VecX() As Double, _
              ByRef VecY() As Double, _
              ByRef a As Double, _
              ByRef b As Double)

    Dim k As Integer, N As Integer
    Dim xm As Double, var As Double, dx As Double
    Dim ym As Double, cov As Double
    
    N = k2 - k1 + 1
    If N < 2 Then Stop   'Données non valides
    xm = 0#: ym = 0#
    For k = k1 To k2
        xm = xm + VecX(k): ym = ym + VecY(k)
    Next k
    xm = xm / N: ym = ym / N
    cov = 0#: var = 0#
    For k = k1 To k2
        dx = VecX(k) - xm
        var = var + dx * dx
        cov = cov + dx * (VecY(k) - ym)
    Next k
    a = cov / var: b = ym - a * xm
End Sub
Function MC_Pol(ByVal k1 As Integer, _
                ByVal k2 As Integer, _
                ByRef VecX() As Double, _
                ByRef VecY() As Double, _
                ByVal m As Integer, _
                ByRef c() As Double) As Integer

    Dim MAT() As Double
    Dim i As Integer, j As Integer, k As Integer
    Dim N As Integer
    Dim u As Double, v As Double, w As Double
    
    If (m < 1) Or (m > k2 - k1) Then Stop
    N = m + 1: ReDim MAT(N, N), c(N)
    For i = 1 To N
        c(i) = 0#: MAT(i, 1) = 0#: MAT(N, i) = 0#
    Next i
    MAT(1, 1) = k2 - k1 + 1
    For k = k1 To k2
        u = VecX(k): v = VecY(k): w = 1#
        c(1) = c(1) + v
        For i = 2 To N
            v = u * v: c(i) = c(i) + v
            w = u * w: MAT(i, 1) = MAT(i, 1) + w
        Next i
        For j = 2 To N
            w = u * w: MAT(N, j) = MAT(N, j) + w
        Next j
    Next k
    For j = 2 To N
        For i = 2 To N
            MAT(i - 1, j) = MAT(i, j - 1)
        Next i
    Next j
    MC_Pol = SYSLIN.Crout(N, MAT(), c())
End Function
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'          Coef. d'approximation Quasi-minimax
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Sub Tchebychev(ByVal nf As Integer, _
               ByVal a As Double, _
               ByVal b As Double, _
               ByVal eps As Double, _
               ByVal N As Integer, _
               ByRef c() As Double, _
               ByRef E() As Double)
               
End Sub
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'          Coef. d'approximation de Fourier
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Sub Fourier(ByVal N As Integer, _
            ByRef y() As Double, _
            ByRef Re() As Double, _
            ByRef Im() As Double)

    Const PI As Double = 3.14159265359
    Dim i As Integer, k As Integer, m As Integer
    Dim p As Double, w  As Double
    Dim c As Double, s  As Double, arg As Double
 
    If (N <= 1) Then Stop
    ReDim Re(1 To N), Im(1 To N)
    m = N / 2
    p = PI / m
    For i = 1 To N
        c = 0#: s = 0#:   w = p * (i - m)
        For k = 1 To N
            arg = w * (k - 1)
            c = c + y(k) * Cos(arg)
            s = s + y(k) * Sin(arg)
        Next k
        Re(i) = c:   Im(i) = s
    Next i
End Sub


