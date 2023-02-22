Attribute VB_Name = "AJUSTE"
Option Explicit
Sub DroiteMC(ByVal k1 As Integer, _
             ByVal k2 As Integer, _
             ByRef x() As Double, _
             ByRef y() As Double, _
             ByRef A As Double, _
             ByRef b As Double)
    
    Dim N As Integer, k As Integer
    Dim xm As Double, ym As Double
    Dim cov As Double, var As Double, e As Double
    
    N = k2 - k1 + 1
    If N < 2 Then Stop
    xm = 0: ym = 0
    For k = k1 To k2
        xm = xm + x(k)
        ym = ym + y(k)
    Next k
    xm = xm / N:  ym = ym / N
    cov = 0#: var = 0#
    For k = k1 To k2
        e = x(k) - xm
        var = var + e * e
        cov = cov + e * (y(k) - ym)
    Next k
    A = cov / var
    b = ym - A * xm
End Sub
Function PolyMc(ByVal k1 As Integer, _
                ByVal k2 As Integer, _
                ByRef x() As Double, _
                ByRef y() As Double, _
                ByVal Deg As Integer, _
                ByRef Coef() As Double) As Integer
                
    Dim N As Integer, k As Integer
    Dim i As Integer, j As Integer
    Dim A() As Double
    
    N = k2 - k1 + 1
    If N < Deg Or Deg < 2 Then Stop
    ReDim A(1 To Deg + 1, 1 To Deg + 1), Coef(1 To Deg + 1)
    For i = 0 To Deg
        For k = k1 To k2
            Coef(i + 1) = Coef(i + 1) + y(k) * (x(k) ^ i)
        Next k
        For j = 0 To Deg
            For k = k1 To k2
                A(i + 1, j + 1) = A(i + 1, j + 1) + x(k) ^ (i + j)
            Next k
        Next j
    Next i
    PolyMc = SYSLIN.LU_FacResol2(Deg + 1, A(), Coef())
End Function

