Attribute VB_Name = "ARCHIVE"
Option Explicit
Sub test_esd()
    
    Dim omar As Integer
    Dim x1 As Double, x2 As Double

    Call TEMP.ESD(0.25, -7.8, x1, x2, omar)
    Debug.Print omar
    Debug.Print x1
    Debug.Print x2

End Sub
Sub test_localise()

    Dim t(1 To 8) As Double
    Dim p As Integer
    
    t(1) = 1.56
    t(2) = 11.7856
    t(3) = -71.56
    t(4) = 11.56
    t(5) = 451.56
    t(6) = 11.56
    t(7) = 31.14756
    t(8) = 21.56
    
    Call TEMP.localise(8, t(), 11.56, 1, p)
    Debug.Print p
    Call TEMP.localise(8, t(), 11.56, -1, p)
    Debug.Print p
    Call TEMP.localise(8, t(), 11.56, 99, p)
    Debug.Print p
    
End Sub
Sub test_Encadre()

    Dim t(1 To 8) As Double
    Dim p1 As Integer
    Dim p2 As Integer
    
    t(1) = -121.56
    t(2) = -111.7856
    t(3) = -71.56
    t(4) = 11.56
    t(5) = 41.56
    t(6) = 411.56
    t(7) = 531.14756
    t(8) = 721.56
    
    Call TEMP.Encadre(8, t(), -230, p1, p2)
    Debug.Print p1
    Debug.Print p2
    
End Sub
Sub test_Remplir()

    Dim x(1 To 21) As Double
    Dim y(1 To 21) As Double
    Dim i As Integer
    
    Call TEMP.Remplir(2, -1, 1, 21, x(), y())
    For i = 1 To 21
        Debug.Print Format(y(i), "0.000E+00   ")
    Next i

End Sub
