Attribute VB_Name = "TEMP"
Option Explicit
Sub ESD(ByVal b As Double, _
        ByVal c As Double, _
        ByRef x1 As Double, _
        ByRef x2 As Double, _
        ByRef sigma As Integer)
        
    'Commentaires
    
    x1 = -0.5 * b
    x2 = x1 * x1 - c
    sigma = Sgn(x2)
    x2 = Sqr(x2 * sigma)
    If sigma > 0 Then
        If x1 > 0 Then x2 = x1 + x2 Else x2 = x1 - x2
        x1 = c / x2
    End If
End Sub
Sub localise(ByVal N As Integer, _
             ByRef t() As Double, _
             ByVal vcible As Double, _
             ByVal SENS As Integer, _
             ByRef p As Integer)
             
    Dim i As Integer, fin As Integer

    If Abs(SENS) <> 1 Or N < 1 Then Stop
    If SENS = 1 Then i = 0 Else i = N + 1
    fin = (N - i) * SENS
    p = 0
    Do
        i = i + SENS
        If t(i) = vcible Then p = i
    Loop Until p <> 0 Or i = fin
End Sub
Sub Encadre(ByVal N As Integer, _
            ByRef x() As Double, _
            ByVal vcible As Double, _
            ByRef p1 As Integer, _
            ByRef p2 As Integer)
    
    Dim m As Integer
    
    If vcible <= x(2) Then
        p1 = 1
    ElseIf vcible >= x(N - 1) Then
        p1 = N - 1
    Else
        p1 = 2
        p2 = N - 1
        Do
            m = (p1 + p2) / 2
            If vcible > x(m) Then p1 = m Else p2 = m
        Loop Until p2 - p1 < 2
    End If
    p2 = p1 + 1
End Sub
Sub Remplir(ByVal nf As Integer, _
            ByVal a As Double, _
            ByVal b As Double, _
            ByVal N As Integer, _
            ByRef x() As Double, _
            ByRef y() As Double)

    Dim h As Double, k As Integer, u As Double

    If N < 2 Then Stop
    h = (b - a) / (N - 1)
    u = a
    For k = 1 To N
        x(k) = u
        y(k) = Fonc(nf, u)
        u = u + h
    Next k
End Sub
